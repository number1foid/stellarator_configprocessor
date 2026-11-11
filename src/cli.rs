//! Command-line interface.
//!
//! A run has two phases. [`prepare`] reads the command line and both input
//! files, allocates every buffer and opens the output files. [`execute`]
//! computes and writes the results and allocates nothing, so the program
//! can lock its allocator between the two.

use crate::consistency::N_TESTS;
use crate::equilibrium::{read_wout, Equilibrium};
use crate::error::{ensure, Code, Error};
use crate::limits::{MAX_ARGS, MAX_FILE_BYTES, MAX_GROUPS, MAX_TEXT};
use crate::makegrid::{read_makegrid, CoilSet, GroupFactors};
use crate::num::{GENERATOR, SCHEMA_VERSION};
use crate::pipeline::{default_settings, new_workspace, run, Outcome, Settings, Workspace};
use crate::report::{failure_text, warning_text, write_json_report, write_text_report};
use crate::schema::{process_file_name, write_config};
use std::fs::{File, OpenOptions};
use std::io::{Seek, Write};
use std::time::Instant;

/// Exit code: success.
pub const EXIT_OK: u8 = 0;
/// Exit code: wrong command line, unreadable input or failed calculation.
pub const EXIT_INPUT: u8 = 2;
/// Exit code: coils and equilibrium are inconsistent.
pub const EXIT_INCONSISTENT: u8 = 3;

const USAGE: &str = "\
usage: stellarator-configprocessor WOUT COILS --wp-radial M --wp-toroidal M
                                   --epseff X --neutron-peakfactor X [options]

Compute the stellarator configuration file (stella_conf.json) that the PROCESS
systems code reads with istell = 6, from a VMEC equilibrium and a MAKEGRID
coil file.

positional arguments:
  WOUT                      VMEC wout file (NetCDF classic format)
  COILS                     MAKEGRID coils file listing all coils

required settings (no defaults):
  --wp-radial M             winding pack width in the coil's in-plane normal
                            direction [m]
  --wp-toroidal M           winding pack width in the toroidal direction [m]
  --epseff X                effective helical ripple (not computed; passed to
                            PROCESS)
  --neutron-peakfactor X    neutron wall load peaking factor (not computed;
                            passed to PROCESS)

output:
  -o, --output FILE         configuration file to write (default:
                            stella_conf.json)
  --process-input IN.DAT    PROCESS input file; the configuration is written
                            next to it under the name PROCESS expects
  --report FILE             report file; .json gives JSON, anything else
                            plain text
  --name NAME               configuration name stored in the file
  -q, --quiet               print nothing on success

reference point:
  --rmajor-ref M            rescale the machine to this plasma major radius [m]
  --bt-ref T                rescale coil currents to this mean field on axis [T]
  --group-current G=FACTOR  multiply the currents of coil group G by FACTOR
                            (VMEC extcur); repeatable

consistency check:
  --match-field             scale all coil currents so the field on axis
                            equals the equilibrium's, instead of checking it
  --allow-inconsistent      continue even if coils and equilibrium appear
                            inconsistent
  --tol-field X             relative tolerance on the field on axis
                            (default 0.03)
  --tol-axis X              tolerance on the field component perpendicular to
                            the axis (default 0.02)
  --tol-linkage X           tolerance on the linked current (default 0.01)

resolution:
  --points N                samples per coil (default 160)
  --wp-filaments NR NT      filaments across the winding pack (default 4 4)

  -h, --help                show this help and exit
  --version                 show the version and exit
";

/// What the command line asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// Compute a configuration.
    Run,
    /// Print the usage text.
    Help,
    /// Print the version.
    Version,
}

/// An option that takes a value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Flag {
    WpRadial,
    WpToroidal,
    Epseff,
    PeakFactor,
    Output,
    ProcessInput,
    Report,
    Name,
    RmajorRef,
    BtRef,
    GroupCurrent,
    TolField,
    TolAxis,
    TolLinkage,
    Points,
    FilamentsRadial,
    FilamentsToroidal,
}

/// The parsed command line.
#[derive(Debug, Clone)]
pub struct Options {
    /// What to do.
    pub action: Action,
    /// Path of the equilibrium file.
    pub wout: String,
    /// Path of the coil file.
    pub coils: String,
    /// Path of the configuration file to write.
    pub output: String,
    /// Path of the report file; empty for none.
    pub report: String,
    /// Configuration name stored in the file.
    pub name: String,
    /// Print nothing on success.
    pub quiet: bool,
    /// Run settings.
    pub settings: Settings,
    /// Factors on the currents of coil groups.
    pub groups: GroupFactors,
}

/// Everything a run needs, ready before the calculation starts.
#[derive(Debug)]
pub struct Session {
    /// The parsed command line.
    pub options: Options,
    /// The equilibrium.
    pub eq: Equilibrium,
    /// The coils.
    pub coil_set: CoilSet,
    /// Every buffer of the calculation.
    pub work: Workspace,
    /// The configuration file, open for writing.
    pub config_file: File,
    /// The report file, open for writing, if one was asked for.
    pub report_file: Option<File>,
}

fn parse_real(text: &str, message: &'static str) -> Result<f64, Error> {
    ensure(!text.is_empty(), Code::Usage, message)?;
    let Ok(value) = text.parse::<f64>() else {
        return Err(Error::new(Code::Usage, message));
    };
    ensure(value.is_finite(), Code::Usage, message)?;
    Ok(value)
}

fn parse_count(text: &str, message: &'static str) -> Result<usize, Error> {
    ensure(!text.is_empty(), Code::Usage, message)?;
    let Ok(value) = text.parse::<usize>() else {
        return Err(Error::new(Code::Usage, message));
    };
    ensure(value <= MAX_FILE_BYTES, Code::Usage, message)?;
    Ok(value)
}

fn parse_group(options: &mut Options, text: &str) -> Result<(), Error> {
    let message = "expected GROUP=FACTOR, for example 1=15000";
    let Some((group, factor)) = text.split_once('=') else {
        return Err(Error::new(Code::Usage, message));
    };
    let index = options.groups.count;
    ensure(index < MAX_GROUPS, Code::Usage, "too many --group-current settings")?;
    let id = parse_real(group, message)?;
    ensure(id.fract() == 0.0, Code::Usage, message)?;
    options.groups.group[index] = id;
    options.groups.factor[index] = parse_real(factor, message)?;
    options.groups.count = index + 1;
    Ok(())
}

/// The value-taking option named `key`, if it is one.
fn flag_for(key: &str) -> Option<Flag> {
    const TABLE: [(&str, Flag); 16] = [
        ("--wp-radial", Flag::WpRadial),
        ("--wp-toroidal", Flag::WpToroidal),
        ("--epseff", Flag::Epseff),
        ("--neutron-peakfactor", Flag::PeakFactor),
        ("--output", Flag::Output),
        ("--process-input", Flag::ProcessInput),
        ("--report", Flag::Report),
        ("--name", Flag::Name),
        ("--rmajor-ref", Flag::RmajorRef),
        ("--bt-ref", Flag::BtRef),
        ("--group-current", Flag::GroupCurrent),
        ("--tol-field", Flag::TolField),
        ("--tol-axis", Flag::TolAxis),
        ("--tol-linkage", Flag::TolLinkage),
        ("--points", Flag::Points),
        ("--wp-filaments", Flag::FilamentsRadial),
    ];
    let mut found = None;
    for (name, flag) in TABLE {
        if name == key {
            found = Some(flag);
        }
    }
    found
}

/// Store the value of a value-taking option.
fn apply_value(options: &mut Options, flag: Flag, text: &str) -> Result<(), Error> {
    ensure(text.len() <= MAX_TEXT, Code::Usage, "an argument is too long")?;
    let number = "an option needs a number";
    let whole = "--wp-filaments needs two whole numbers";
    let settings = &mut options.settings;
    match flag {
        Flag::WpRadial => settings.wp_radial = parse_real(text, number)?,
        Flag::WpToroidal => settings.wp_toroidal = parse_real(text, number)?,
        Flag::Epseff => settings.epseff = parse_real(text, number)?,
        Flag::PeakFactor => settings.neutron_peakfactor = parse_real(text, number)?,
        Flag::RmajorRef => settings.rmajor_ref = Some(parse_real(text, number)?),
        Flag::BtRef => settings.bt_ref = Some(parse_real(text, number)?),
        Flag::TolField => settings.tolerances.field = parse_real(text, number)?,
        Flag::TolAxis => settings.tolerances.axis_misalignment = parse_real(text, number)?,
        Flag::TolLinkage => settings.tolerances.linkage = parse_real(text, number)?,
        Flag::Points => settings.n_points = parse_count(text, "--points needs a whole number")?,
        Flag::FilamentsRadial => settings.wp_filaments.0 = parse_count(text, whole)?,
        Flag::FilamentsToroidal => settings.wp_filaments.1 = parse_count(text, whole)?,
        Flag::Output => {
            ensure(
                options.output.is_empty(),
                Code::Usage,
                "give either --output or --process-input, once",
            )?;
            options.output.push_str(text);
        }
        Flag::ProcessInput => {
            ensure(
                options.output.is_empty(),
                Code::Usage,
                "give either --output or --process-input, once",
            )?;
            options.output = process_file_name(text)?;
        }
        Flag::Report => options.report.push_str(text),
        Flag::Name => options.name.push_str(text),
        Flag::GroupCurrent => parse_group(options, text)?,
    }
    Ok(())
}

/// Handle an argument that is not the value of an option. Returns the
/// option that now waits for its value, if any.
fn apply_argument(
    options: &mut Options,
    argument: &str,
    positional: &mut usize,
) -> Result<Option<Flag>, Error> {
    ensure(argument.len() <= MAX_TEXT, Code::Usage, "an argument is too long")?;
    let (head, tail) = argument.split_once('=').unwrap_or((argument, ""));
    let key = if head == "-o" { "--output" } else { head };
    let mut waiting = None;
    if let Some(flag) = flag_for(key) {
        if argument.contains('=') {
            apply_value(options, flag, tail)?;
        } else {
            waiting = Some(flag);
        }
    } else if argument == "-h" || argument == "--help" {
        options.action = Action::Help;
    } else if argument == "--version" {
        options.action = Action::Version;
    } else if argument == "-q" || argument == "--quiet" {
        options.quiet = true;
    } else if argument == "--match-field" {
        options.settings.match_field = true;
    } else if argument == "--allow-inconsistent" {
        options.settings.allow_inconsistent = true;
    } else {
        ensure(!argument.starts_with('-'), Code::Usage, "unknown option; see --help")?;
        ensure(*positional < 2, Code::Usage, "too many arguments; expected WOUT and COILS")?;
        if *positional == 0 {
            options.wout.push_str(argument);
        } else {
            options.coils.push_str(argument);
        }
        *positional += 1;
    }
    Ok(waiting)
}

/// Check that the command line is complete and fill in the defaults.
fn finish_options(options: &mut Options, positional: usize) -> Result<(), Error> {
    if options.action != Action::Run {
        return Ok(());
    }
    ensure(positional == 2, Code::Usage, "expected two arguments: WOUT and COILS; see --help")?;
    let settings = options.settings;
    ensure(
        !(settings.wp_radial.is_nan() || settings.wp_toroidal.is_nan()),
        Code::Usage,
        "--wp-radial and --wp-toroidal are required",
    )?;
    ensure(
        !(settings.epseff.is_nan() || settings.neutron_peakfactor.is_nan()),
        Code::Usage,
        "--epseff and --neutron-peakfactor are required",
    )?;
    if options.output.is_empty() {
        options.output.push_str("stella_conf.json");
    }
    if options.name.is_empty() {
        options.name.push_str("stellarator_configprocessor");
    }
    Ok(())
}

/// Parse the command-line arguments (without the program name).
///
/// # Errors
/// Returns an error of class [`Code::Usage`] if an option is unknown, a
/// value is missing or malformed, or a required setting is absent.
pub fn parse_args<I: Iterator<Item = String>>(arguments: I) -> Result<Options, Error> {
    let mut options = Options {
        action: Action::Run,
        wout: String::new(),
        coils: String::new(),
        output: String::new(),
        report: String::new(),
        name: String::new(),
        quiet: false,
        settings: default_settings((f64::NAN, f64::NAN), f64::NAN, f64::NAN),
        groups: GroupFactors { count: 0, group: [0.0; MAX_GROUPS], factor: [1.0; MAX_GROUPS] },
    };
    let mut waiting: Option<Flag> = None;
    let (mut positional, mut seen) = (0, 0);
    for argument in arguments.take(MAX_ARGS + 1) {
        seen += 1;
        if let Some(flag) = waiting {
            apply_value(&mut options, flag, &argument)?;
            waiting =
                if flag == Flag::FilamentsRadial { Some(Flag::FilamentsToroidal) } else { None };
        } else {
            waiting = apply_argument(&mut options, &argument, &mut positional)?;
        }
    }
    ensure(seen <= MAX_ARGS, Code::Usage, "too many arguments")?;
    ensure(waiting.is_none(), Code::Usage, "the last option needs a value")?;
    finish_options(&mut options, positional)?;
    Ok(options)
}

/// Open a file for writing without discarding what it holds, so that a
/// failed run leaves an existing file untouched.
fn open_output(path: &str) -> Result<File, Error> {
    ensure(!path.is_empty(), Code::Usage, "an output file needs a name")?;
    ensure(path.len() <= MAX_TEXT, Code::Usage, "an output file name is too long")?;
    let Ok(file) = OpenOptions::new().write(true).create(true).truncate(false).open(path) else {
        return Err(Error::new(Code::Io, "cannot open an output file for writing"));
    };
    Ok(file)
}

/// Read the command line and the input files, allocate every buffer and
/// open the output files. Returns `None` if only help or the version was
/// asked for, after writing it to `out`.
///
/// # Errors
/// Returns an error if the command line is wrong, an input file cannot be
/// read or is unusable, or an output file cannot be opened.
pub fn prepare<I: Iterator<Item = String>, W: Write>(
    arguments: I,
    out: &mut W,
) -> Result<Option<Session>, Error> {
    let options = parse_args(arguments)?;
    if options.action == Action::Help {
        write!(out, "{USAGE}\nWrites schema {SCHEMA_VERSION}. See docs/user_guide.md.\n")?;
        return Ok(None);
    }
    if options.action == Action::Version {
        writeln!(out, "{GENERATOR}")?;
        return Ok(None);
    }
    let Ok(bytes) = std::fs::read(&options.wout) else {
        return Err(Error::new(Code::Equilibrium, "cannot read the equilibrium file"));
    };
    let Ok(text) = std::fs::read_to_string(&options.coils) else {
        return Err(Error::new(Code::CoilFile, "cannot read the coil file"));
    };
    ensure(bytes.len() <= MAX_FILE_BYTES, Code::Limit, "the equilibrium file is too large")?;
    let eq = read_wout(&bytes)?;
    let coil_set = read_makegrid(&text, &options.groups)?;
    let work = new_workspace(&eq, &coil_set, &options.settings)?;
    let matching = work.coils.n_coils == coil_set.currents.len();
    ensure(matching, Code::Numeric, "the workspace does not match the coil set")?;
    let config_file = open_output(&options.output)?;
    let report_file =
        if options.report.is_empty() { None } else { Some(open_output(&options.report)?) };
    Ok(Some(Session { options, eq, coil_set, work, config_file, report_file }))
}

/// Whether a path ends in `.json`, in any letter case.
fn is_json(path: &str) -> bool {
    let bytes = path.as_bytes();
    bytes.len() >= 5 && bytes[bytes.len() - 5..].eq_ignore_ascii_case(b".json")
}

/// Write the results of a successful run.
fn write_results<W: Write>(
    session: &mut Session,
    outcome: &Outcome,
    runtime: f64,
    out: &mut W,
) -> Result<(), Error> {
    ensure(runtime.is_finite(), Code::Numeric, "the run time is not valid")?;
    let options = &session.options;
    let inputs = (options.wout.as_str(), options.coils.as_str(), options.name.as_str());
    write_config(&mut session.config_file, &outcome.config, &options.name)?;
    let end = session.config_file.stream_position()?;
    session.config_file.set_len(end)?;
    ensure(end > 0, Code::Io, "the configuration file is empty")?;
    if !options.quiet {
        writeln!(out)?;
        write_text_report(out, inputs, &options.settings, outcome, &session.work, runtime)?;
        writeln!(out, "wrote {}", options.output)?;
    }
    if let Some(file) = session.report_file.as_mut() {
        if is_json(&options.report) {
            write_json_report(file, inputs, &options.settings, outcome, &session.work, runtime)?;
        } else {
            write_text_report(file, inputs, &options.settings, outcome, &session.work, runtime)?;
        }
        let end = file.stream_position()?;
        file.set_len(end)?;
        if !options.quiet {
            writeln!(out, "wrote {}", options.report)?;
        }
    }
    Ok(())
}

/// Compute and write the results; print warnings to `err`.
fn compute<W: Write, E: Write>(
    session: &mut Session,
    out: &mut W,
    err: &mut E,
) -> Result<(), Error> {
    let start = Instant::now();
    let settings = session.options.settings;
    ensure(session.options.action == Action::Run, Code::Usage, "nothing to compute")?;
    let outcome = if session.options.quiet {
        run(
            &mut session.eq,
            &mut session.coil_set,
            &settings,
            &mut session.work,
            &mut std::io::sink(),
        )?
    } else {
        run(&mut session.eq, &mut session.coil_set, &settings, &mut session.work, out)?
    };
    let runtime = start.elapsed().as_secs_f64();
    write_results(session, &outcome, runtime, out)?;
    for (index, active) in outcome.warnings.iter().enumerate() {
        if *active {
            write!(err, "warning: ")?;
            warning_text(err, &settings, &outcome, index)?;
            writeln!(err)?;
        }
    }
    ensure(outcome.config.assigned[1], Code::Schema, "the configuration is incomplete")?;
    Ok(())
}

/// Write an error, with the failed consistency tests if that is the cause,
/// to `err`. Returns the exit code.
///
/// # Errors
/// Returns an error if writing fails.
pub fn report_failure<E: Write>(
    err: &mut E,
    error: &Error,
    session: Option<&Session>,
) -> Result<u8, Error> {
    ensure(!error.message.is_empty(), Code::Numeric, "an error has no message")?;
    writeln!(err, "error: {error}")?;
    if error.code != Code::Inconsistent {
        return Ok(EXIT_INPUT);
    }
    if let Some(Session { options, work, .. }) = session {
        if let Some(check) = work.last_check {
            for index in 0..N_TESTS {
                if check.failed[index] {
                    write!(err, "  - ")?;
                    failure_text(err, &check, &options.settings.tolerances, index)?;
                    writeln!(err)?;
                }
            }
        }
    }
    writeln!(
        err,
        "Use --allow-inconsistent to continue anyway, or --match-field if only the current \
         scale differs."
    )?;
    Ok(EXIT_INCONSISTENT)
}

/// Run the calculation of a prepared session and write its results.
/// Returns the exit code of the program. This function allocates nothing.
pub fn execute<W: Write, E: Write>(session: &mut Session, out: &mut W, err: &mut E) -> u8 {
    match compute(session, out, err) {
        Ok(()) => EXIT_OK,
        Err(error) => report_failure(err, &error, Some(session)).unwrap_or(EXIT_INPUT),
    }
}
