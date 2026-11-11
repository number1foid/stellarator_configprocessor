//! Coil file reader, schema writer and command line.

use crate::helpers::{close, failure, NO_FACTORS};
use stellarator_configprocessor::cli::{parse_args, Action};
use stellarator_configprocessor::error::{Code, Error};
use stellarator_configprocessor::makegrid::{read_makegrid, scale_coils, GroupFactors};
use stellarator_configprocessor::schema::{
    process_file_name, set_value, validate, write_config, write_json_text, Kind, EMPTY_CONFIG,
    KEYS, N_KEYS,
};

const SQUARE: &str = "\
&nml  ignored = 1 /
periods 3
begin filament
mirror NIL
 0.0 0.0 0.0  1.0D+03
 1.0 0.0 0.0  1.0D+03
 1.0 1.0 0.0  1.0D+03
 0.0 1.0 0.0  1.0D+03
 0.0 0.0 0.0  0.0  2 square
 3.0 0.0 0.5  -2.5E2  7 ring
end
";

fn arguments(words: &[&str]) -> std::vec::IntoIter<String> {
    let owned: Vec<String> = words.iter().map(|word| (*word).to_owned()).collect();
    owned.into_iter()
}

#[test]
fn square_and_ring_are_read() -> Result<(), Error> {
    let set = read_makegrid(SQUARE, &NO_FACTORS)?;
    assert_eq!(set.periods, 3);
    assert_eq!(set.offsets, vec![0, 4, 4 + 128]);
    close(set.currents[0], 1000.0, 0.0, "current in Fortran notation");
    close(set.currents[1], -250.0, 0.0, "current of the one-line coil");
    close(set.points[4][0], 3.0, 1.0e-15, "first point of the circle");
    close(set.points[4 + 32][1], 3.0, 1.0e-12, "quarter point of the circle");
    close(set.points[4 + 32][2], 0.5, 0.0, "height of the circle");
    let notes = set.notes;
    assert_eq!((notes.zero_current, notes.not_closed, notes.varying_current), (0, 0, 0));
    Ok(())
}

#[test]
fn group_factors_and_scaling_apply() -> Result<(), Error> {
    let mut factors = GroupFactors { count: 2, group: [0.0; 32], factor: [1.0; 32] };
    factors.group[0] = 2.0;
    factors.factor[0] = 3.0;
    factors.group[1] = 7.0;
    factors.factor[1] = 0.0;
    let mut set = read_makegrid(SQUARE, &factors)?;
    assert_eq!(set.currents.len(), 1, "the ring has zero current and is left out");
    close(set.currents[0], 3000.0, 0.0, "scaled group current");
    scale_coils(&mut set, 2.0, -0.5)?;
    close(set.points[2][1], 2.0, 0.0, "scaled coordinate");
    close(set.currents[0], -1500.0, 0.0, "scaled current");
    assert!(scale_coils(&mut set, 0.0, 1.0).is_err());
    Ok(())
}

#[test]
fn irregular_coils_are_noted() -> Result<(), Error> {
    let open =
        "periods 1\nbegin filament\n0 0 0 5\n1 0 0 5\n1 1 0 6\n0 1 0 5\n0 0.2 0 1 1 open\nend\n";
    let set = read_makegrid(open, &NO_FACTORS)?;
    assert_eq!(set.offsets, vec![0, 5], "the closing point is kept for an open coil");
    assert_eq!(set.notes.not_closed, 1);
    assert_eq!(set.notes.varying_current, 1);
    assert_eq!(set.notes.closing_current, 1);
    let idle = "periods 1\n0 0 0 0\n1 0 0 0\n1 1 0 0\n0 0 0 0 1 idle\n1 0 0 5 1 ring\nend\n";
    let set = read_makegrid(idle, &NO_FACTORS)?;
    assert_eq!((set.currents.len(), set.notes.zero_current), (1, 1));
    Ok(())
}

#[test]
fn malformed_coil_files_are_rejected() {
    let cases = [
        ("begin filament\n0 0 0 1\nend\n", "no 'periods' line", 0),
        ("periods two\n", "cannot parse numbers", 1),
        ("periods 2.5\n", "bad 'periods' line", 1),
        ("periods 2\nbegin filament\n1 0 0 5\n2 0 0 5\nend\n", "ends inside a coil", 0),
        ("periods 2\nbegin filament\n1 0 x 5\nend\n", "cannot parse numbers", 3),
        ("periods 2\n1 0 0\nend\n", "expected 'x y z current'", 2),
        ("periods 2\n1 0 0 5\n1 0 0 5\n1 0 0 0 1 dot\nend\n", "no extent", 4),
        ("periods 2\n0 0 0 5\n1 0 0 5\n0 0 0 0 1 line\nend\n", "at least 3 distinct points", 4),
        ("periods 2\n0 0 0 5\n1 0 0 5\n1 1 0 5\n0 0 0 0 one sq\nend\n", "bad group id", 5),
        ("periods 2\nbegin filament\nend\n", "no coils found", 0),
        (
            "periods 2\n0 0 0 0\n1 0 0 0\n1 1 0 0\n0 0 0 0 1 idle\nend\n",
            "no coil carries current",
            0,
        ),
        ("periods 2\n0 0 0 inf\nend\n", "non-finite number", 2),
    ];
    for (text, message, line) in cases {
        let error = failure(read_makegrid(text, &NO_FACTORS), message);
        assert_eq!(error.code, Code::CoilFile, "{message}");
        assert!(error.message.contains(message), "expected `{message}`, got `{}`", error.message);
        assert_eq!(error.line, line, "{message}");
    }
}

#[test]
fn schema_has_the_keys_process_reads() {
    assert_eq!(KEYS.len(), N_KEYS);
    assert_eq!(KEYS[0].name, "name");
    assert_eq!(KEYS[0].kind, Kind::Text);
    assert_eq!(KEYS[N_KEYS - 1].name, "neutron_peakfactor");
    for (index, key) in KEYS.iter().enumerate() {
        let repeated = KEYS[..index].iter().any(|other| other.name == key.name);
        assert!(!repeated, "key `{}` appears twice", key.name);
        assert!(!key.description.is_empty() && !key.method.is_empty());
    }
}

#[test]
fn configuration_is_validated_and_written() -> Result<(), Error> {
    let mut config = EMPTY_CONFIG;
    assert_eq!(failure(validate(&config), "empty").code, Code::Schema);
    for key in &KEYS {
        if key.kind != Kind::Text {
            set_value(&mut config, key.name, if key.kind == Kind::Count { 5.0 } else { 1.25 })?;
        }
    }
    validate(&config)?;
    assert!(set_value(&mut config, "bt_ref", 2.0).is_err(), "a key may be set only once");
    assert!(set_value(&mut config, "no_such_key", 2.0).is_err());
    assert!(set_value(&mut config, "name", 2.0).is_err());
    let mut text = Vec::new();
    write_config(&mut text, &config, "a \"quoted\" name")?;
    let text = String::from_utf8_lossy(&text);
    assert!(text.starts_with("{\n  \"name\": \"a \\\"quoted\\\" name\",\n  \"symmetry\": 5,\n"));
    assert!(text.contains("\"rmajor_ref\": 1.25,"));
    assert!(text.ends_with("\"_generator\": \"stellarator_configprocessor 0.2.0\"\n}\n"));
    assert_eq!(text.matches(": ").count(), N_KEYS + 2);

    let mut broken = config;
    broken.values[3] = -1.0;
    assert!(validate(&broken).is_err(), "rmajor_ref must be positive");
    broken.values[3] = f64::NAN;
    assert!(validate(&broken).is_err(), "values must be finite");
    broken = config;
    broken.values[1] = 2.5;
    assert!(validate(&broken).is_err(), "symmetry must be a whole number");
    Ok(())
}

#[test]
fn json_text_is_escaped() -> Result<(), Error> {
    let mut text = Vec::new();
    write_json_text(&mut text, "C:\\data\\wout \"x\"\tend")?;
    assert_eq!(String::from_utf8_lossy(&text), "\"C:\\\\data\\\\wout \\\"x\\\"\\u0009end\"");
    Ok(())
}

#[test]
fn process_file_names() -> Result<(), Error> {
    assert_eq!(process_file_name("run1.IN.DAT")?, "run1.stella_conf.json");
    assert_eq!(process_file_name("out/IN.DAT")?, "out/stella_conf.json");
    assert_eq!(failure(process_file_name("run1.txt"), "name").code, Code::Usage);
    Ok(())
}

#[test]
fn command_line_is_parsed() -> Result<(), Error> {
    let options = parse_args(arguments(&[
        "wout.nc",
        "coils.x",
        "--wp-radial",
        "0.6",
        "--wp-toroidal=0.5",
        "--epseff",
        "0.01",
        "--neutron-peakfactor",
        "1.6",
        "-o",
        "out.json",
        "--bt-ref",
        "5.5",
        "--points",
        "200",
        "--wp-filaments",
        "3",
        "5",
        "--group-current",
        "2=1.5",
        "--match-field",
        "-q",
    ]))?;
    assert_eq!(options.action, Action::Run);
    assert_eq!((options.wout.as_str(), options.coils.as_str()), ("wout.nc", "coils.x"));
    assert_eq!(options.output, "out.json");
    assert_eq!(options.name, "stellarator_configprocessor");
    let settings = options.settings;
    close(settings.wp_radial, 0.6, 0.0, "wp_radial");
    close(settings.wp_toroidal, 0.5, 0.0, "wp_toroidal");
    close(settings.bt_ref.unwrap_or(0.0), 5.5, 0.0, "bt_ref");
    assert!(settings.rmajor_ref.is_none() && settings.match_field && options.quiet);
    assert_eq!((settings.n_points, settings.wp_filaments), (200, (3, 5)));
    assert_eq!(options.groups.count, 1);
    close(options.groups.factor[0], 1.5, 0.0, "group factor");

    let named = parse_args(arguments(&[
        "a",
        "b",
        "--wp-radial=1",
        "--wp-toroidal=1",
        "--epseff=0",
        "--neutron-peakfactor=1",
        "--process-input",
        "runs/big.IN.DAT",
    ]))?;
    assert_eq!(named.output, "runs/big.stella_conf.json");
    assert_eq!(parse_args(arguments(&["--help"]))?.action, Action::Help);
    assert_eq!(parse_args(arguments(&["--version"]))?.action, Action::Version);
    Ok(())
}

#[test]
fn wrong_command_lines_are_rejected() {
    let complete =
        ["a", "b", "--wp-radial=1", "--wp-toroidal=1", "--epseff=0", "--neutron-peakfactor=1"];
    let cases: [&[&str]; 8] = [
        &["a", "b", "--wp-radial=1", "--wp-toroidal=1", "--epseff=0"],
        &["a", "--wp-radial=1", "--wp-toroidal=1", "--epseff=0", "--neutron-peakfactor=1"],
        &["a", "b", "c", "--wp-radial=1"],
        &["a", "b", "--no-such-option"],
        &["a", "b", "--wp-radial", "wide"],
        &["a", "b", "--points"],
        &["a", "b", "--group-current", "two"],
        &["a", "b", "-o", "x.json", "--process-input", "y.IN.DAT"],
    ];
    for words in cases {
        let error = failure(parse_args(arguments(words)), "command line");
        assert_eq!(error.code, Code::Usage, "{words:?}: {}", error.message);
    }
    assert!(parse_args(arguments(&complete)).is_ok());
}
