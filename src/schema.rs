//! The `stella_conf.json` schema read by PROCESS, and its writer.
//!
//! See `docs/schema.md` for the derivation of every entry from the PROCESS
//! source. Key spelling follows the example file shipped with PROCESS;
//! PROCESS itself lower-cases keys on reading.

use crate::error::{ensure, Code, Error};
use crate::limits::MAX_TEXT;
use crate::num::{GENERATOR, SCHEMA_VERSION};
use std::io::Write;

const BAD: Code = Code::Schema;

/// Kind of value stored under a key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A text.
    Text,
    /// A whole number.
    Count,
    /// A real number.
    Real,
}

/// One entry of the configuration file.
#[derive(Debug, Clone, Copy)]
pub struct Key {
    /// Name in the file.
    pub name: &'static str,
    /// Kind of value.
    pub kind: Kind,
    /// Unit, empty if dimensionless.
    pub unit: &'static str,
    /// What the value is.
    pub description: &'static str,
    /// How the value is obtained.
    pub method: &'static str,
    /// Whether the value must be greater than zero.
    pub positive: bool,
}

const fn key(
    name: &'static str,
    kind: Kind,
    unit: &'static str,
    description: &'static str,
    method: &'static str,
    positive: bool,
) -> Key {
    Key { name, kind, unit, description, method, positive }
}

/// Number of entries of the configuration file.
pub const N_KEYS: usize = 35;

/// Every entry PROCESS reads, in file order.
pub const KEYS: [Key; N_KEYS] = [
    key("name", Kind::Text, "", "Configuration name", "run setting", false),
    key("symmetry", Kind::Count, "", "Number of field periods", "from the equilibrium", true),
    key(
        "coilspermodule",
        Kind::Count,
        "",
        "Coils per field period",
        "coil count divided by the number of field periods",
        true,
    ),
    key(
        "rmajor_ref",
        Kind::Real,
        "m",
        "Reference plasma major radius",
        "volume / (2 pi x mean cross-section) of the last closed flux surface",
        true,
    ),
    key(
        "rminor_ref",
        Kind::Real,
        "m",
        "Reference plasma minor radius",
        "radius of a circle with the mean cross-section",
        true,
    ),
    key("aspect_ref", Kind::Real, "", "Reference aspect ratio", "rmajor_ref / rminor_ref", true),
    key(
        "bt_ref",
        Kind::Real,
        "T",
        "Reference field on axis",
        "arc-length mean of |B| of the coils along the magnetic axis",
        true,
    ),
    key(
        "vol_plasma",
        Kind::Real,
        "m^3",
        "Plasma volume",
        "integral over the last closed flux surface (Lion 2021, eq. 3)",
        true,
    ),
    key(
        "plasma_surface",
        Kind::Real,
        "m^2",
        "Plasma surface area",
        "integral over the last closed flux surface (Lion 2021, eq. 4)",
        true,
    ),
    key(
        "coil_rmajor",
        Kind::Real,
        "m",
        "Mean coil major radius",
        "arc-length-weighted mean cylindrical radius of the filaments",
        true,
    ),
    key(
        "coil_rminor",
        Kind::Real,
        "m",
        "Mean coil minor radius",
        "arc-length-weighted mean distance from the circle R = coil_rmajor",
        true,
    ),
    key(
        "coillength",
        Kind::Real,
        "m",
        "Total length of all coils",
        "sum of filament lengths",
        true,
    ),
    key(
        "coilsurface",
        Kind::Real,
        "m^2",
        "Area of the toroidal surface carrying the coils",
        "mean surface of revolution of the filaments: 2 pi x coil_rmajor x mean coil length",
        true,
    ),
    key(
        "maximal_coil_height",
        Kind::Real,
        "m",
        "Full height of the coil set",
        "max(z) - min(z) of the filaments",
        true,
    ),
    key(
        "dmin",
        Kind::Real,
        "m",
        "Minimum distance between two coil filaments",
        "shortest distance between any two different filaments",
        true,
    ),
    key(
        "max_portsize_width",
        Kind::Real,
        "m",
        "Toroidal port width",
        "largest distance between the centroids of neighbouring coils",
        true,
    ),
    key(
        "min_bend_radius",
        Kind::Real,
        "m",
        "Minimum bending radius",
        "1 / max curvature of a periodic spline through each filament (Lion 2021, eq. 47)",
        true,
    ),
    key(
        "min_plasma_coil_distance",
        Kind::Real,
        "m",
        "Minimum plasma-coil distance",
        "shortest distance between the last closed flux surface and any filament",
        true,
    ),
    key(
        "derivative_min_LCFS_coils_dist",
        Kind::Real,
        "",
        "Derivative of the plasma-coil distance with respect to the minor radius",
        "linear fit over flux surfaces near the boundary (Lion 2021, eq. 58)",
        false,
    ),
    key(
        "WP_area",
        Kind::Real,
        "m^2",
        "Reference winding pack area",
        "radial x toroidal width (run settings)",
        true,
    ),
    key("WP_ratio", Kind::Real, "", "Winding pack radial / toroidal width", "run settings", true),
    key(
        "i0",
        Kind::Real,
        "MA",
        "Current per coil for bt_ref at rmajor_ref",
        "mean absolute coil current",
        true,
    ),
    key(
        "WP_bmax",
        Kind::Real,
        "T",
        "Peak field in the reference winding pack",
        "maximum of |B| over the winding pack of all coils (Lion 2021, eq. 38)",
        true,
    ),
    key(
        "a1",
        Kind::Real,
        "",
        "Peak-field fit parameter",
        "least-squares fit of the PROCESS formula over several winding pack sizes",
        false,
    ),
    key(
        "a2",
        Kind::Real,
        "",
        "Peak-field fit parameter (winding pack size term)",
        "least-squares fit of the PROCESS formula over several winding pack sizes",
        false,
    ),
    key(
        "inductance",
        Kind::Real,
        "H",
        "Effective inductance of the coil set",
        "2 x magnetic energy / i0^2",
        true,
    ),
    key(
        "max_force_density",
        Kind::Real,
        "MN/m^3",
        "Maximum mean force density in the winding pack",
        "max_force_density_MNm / WP_area",
        true,
    ),
    key(
        "max_force_density_MNm",
        Kind::Real,
        "MN/m",
        "Maximum force per unit length",
        "maximum of |F'| over all coils (Lion 2021, eq. 40)",
        true,
    ),
    key(
        "max_lateral_force_density",
        Kind::Real,
        "MN/m^3",
        "Maximum lateral force density",
        "maximum toroidal component of F' / WP_area",
        false,
    ),
    key(
        "max_radial_force_density",
        Kind::Real,
        "MN/m^3",
        "Maximum radial force density",
        "maximum in-plane normal component of F' / WP_area",
        true,
    ),
    key(
        "centering_force_max_MN",
        Kind::Real,
        "MN",
        "Largest centering force on a coil",
        "net coil force towards the torus axis",
        false,
    ),
    key(
        "centering_force_min_MN",
        Kind::Real,
        "MN",
        "Smallest centering force on a coil",
        "net coil force towards the torus axis; negative is outwards",
        false,
    ),
    key(
        "centering_force_avg_MN",
        Kind::Real,
        "MN",
        "Mean centering force per coil",
        "average over all coils",
        false,
    ),
    key("epseff", Kind::Real, "", "Effective helical ripple", "run setting, not computed", false),
    key(
        "neutron_peakfactor",
        Kind::Real,
        "",
        "Neutron wall load peaking factor",
        "run setting, not computed",
        true,
    ),
];

/// The numbers of one configuration, in the order of [`KEYS`]. The entry of
/// a text key is unused.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Value of each key.
    pub values: [f64; N_KEYS],
    /// Whether each key has been given a value.
    pub assigned: [bool; N_KEYS],
}

/// A configuration without any value.
pub const EMPTY_CONFIG: Config = Config { values: [0.0; N_KEYS], assigned: [false; N_KEYS] };

/// Store the number `value` under the key `name`.
///
/// # Errors
/// Returns an error if the schema has no such key, the key holds a text, or
/// the key already has a value.
pub fn set_value(config: &mut Config, name: &str, value: f64) -> Result<(), Error> {
    let mut found = N_KEYS;
    for (index, entry) in KEYS.iter().enumerate() {
        if entry.name == name {
            found = index;
        }
    }
    ensure(found < N_KEYS, BAD, "unknown configuration key")?;
    ensure(KEYS[found].kind != Kind::Text, BAD, "a text key cannot hold a number")?;
    ensure(!config.assigned[found], BAD, "a configuration key was set twice")?;
    config.values[found] = value;
    config.assigned[found] = true;
    Ok(())
}

/// Check that every key has a valid value.
///
/// # Errors
/// Returns an error if a key is missing, a value is not finite, a whole
/// number has a fraction, or a value that must be positive is not.
pub fn validate(config: &Config) -> Result<(), Error> {
    for (index, entry) in KEYS.iter().enumerate() {
        if entry.kind != Kind::Text {
            let value = config.values[index];
            ensure(config.assigned[index], BAD, "invalid configuration: a key is missing")?;
            ensure(value.is_finite(), BAD, "invalid configuration: a value is not finite")?;
            ensure(
                entry.kind != Kind::Count || value.fract() == 0.0,
                BAD,
                "invalid configuration: a count is not a whole number",
            )?;
            ensure(
                !entry.positive || value > 0.0,
                BAD,
                "invalid configuration: a value that must be positive is not",
            )?;
        }
    }
    Ok(())
}

/// Write `text` as a JSON string, with quotes and escapes.
///
/// # Errors
/// Returns an error if the text is too long or writing fails.
pub fn write_json_text<W: Write>(out: &mut W, text: &str) -> Result<(), Error> {
    ensure(!text.contains('\u{0}'), BAD, "a text contains a NUL character")?;
    ensure(text.len() <= MAX_TEXT, Code::Limit, "a text is too long to be written")?;
    write!(out, "\"")?;
    for character in text.chars().take(MAX_TEXT) {
        if character == '"' || character == '\\' {
            write!(out, "\\{character}")?;
        } else if character < ' ' {
            write!(out, "\\u{:04x}", u32::from(character))?;
        } else {
            write!(out, "{character}")?;
        }
    }
    write!(out, "\"")?;
    Ok(())
}

/// Write `value` as a JSON number; `null` if it is not finite.
///
/// # Errors
/// Returns an error if writing fails.
pub fn write_json_real<W: Write>(out: &mut W, value: f64) -> Result<(), Error> {
    if value.is_finite() {
        write!(out, "{value:?}")?;
    } else {
        write!(out, "null")?;
    }
    Ok(())
}

/// Validate the configuration and write the file PROCESS reads with
/// `istell = 6`.
///
/// # Errors
/// Returns an error if the configuration is invalid or writing fails.
pub fn write_config<W: Write>(out: &mut W, config: &Config, name: &str) -> Result<(), Error> {
    ensure(name.len() <= MAX_TEXT, Code::Limit, "the configuration name is too long")?;
    validate(config)?;
    ensure(!name.is_empty(), BAD, "the configuration needs a name")?;
    writeln!(out, "{{")?;
    for (index, entry) in KEYS.iter().enumerate() {
        write!(out, "  \"{}\": ", entry.name)?;
        if entry.kind == Kind::Text {
            write_json_text(out, name)?;
        } else if entry.kind == Kind::Count {
            write!(out, "{:.0}", config.values[index])?;
        } else {
            write_json_real(out, config.values[index])?;
        }
        writeln!(out, ",")?;
    }
    writeln!(out, "  \"_schema_version\": \"{SCHEMA_VERSION}\",")?;
    writeln!(out, "  \"_generator\": \"{GENERATOR}\"")?;
    writeln!(out, "}}")?;
    Ok(())
}

/// Name PROCESS expects for the configuration of a given input file.
///
/// PROCESS opens `<prefix>stella_conf.json` where the prefix is the input
/// file name without `IN.DAT`: `run1.IN.DAT` -> `run1.stella_conf.json`.
///
/// # Errors
/// Returns an error if the input file name does not end in `IN.DAT`.
pub fn process_file_name(input_file: &str) -> Result<String, Error> {
    ensure(input_file.len() <= MAX_TEXT, Code::Limit, "the PROCESS input file name is too long")?;
    ensure(
        input_file.ends_with("IN.DAT"),
        Code::Usage,
        "a PROCESS input file name ends in 'IN.DAT'",
    )?;
    let mut name = String::from(&input_file[..input_file.len() - "IN.DAT".len()]);
    name.push_str("stella_conf.json");
    Ok(name)
}
