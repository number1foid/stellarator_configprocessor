//! Reader for the `NetCDF` classic file format (CDF-1 and CDF-2), which VMEC
//! uses for its `wout` files.
//!
//! Only what the equilibrium reader needs is implemented: locating a
//! fixed-size variable by name and reading its numbers. Files in the
//! HDF5-based `NetCDF-4` format are rejected with a clear message.

use crate::error::{ensure, Code, Error};
use crate::limits::{MAX_FILE_BYTES, MAX_NC_DIMS, MAX_NC_ENTRIES, MAX_RANK};

const TAG_DIMENSION: u32 = 10;
const TAG_VARIABLE: u32 = 11;
const TAG_ATTRIBUTE: u32 = 12;
const TYPE_INT: u32 = 4;
const TYPE_FLOAT: u32 = 5;
const TYPE_DOUBLE: u32 = 6;
const BAD: Code = Code::Equilibrium;

/// The parts of the file header needed to locate variables.
#[derive(Debug, Clone, Copy)]
pub struct Header {
    /// Whether data offsets are stored in 64 bits (CDF-2).
    wide_offsets: bool,
    /// Number of dimensions declared in the file.
    n_dims: usize,
    /// Length of each dimension; zero marks the record dimension.
    dim_len: [usize; MAX_NC_DIMS],
    /// Byte position of the first variable entry.
    first_variable: usize,
    /// Number of variables declared in the file.
    n_variables: usize,
}

/// Where the numbers of one variable are stored.
#[derive(Debug, Clone, Copy)]
pub struct Variable {
    /// Number type code of the file format.
    pub kind: u32,
    /// Byte position of the first value.
    pub begin: usize,
    /// Number of values.
    pub count: usize,
    /// Number of dimensions.
    pub rank: usize,
    /// Length of each dimension.
    pub shape: [usize; MAX_RANK],
}

/// A read position in the file.
#[derive(Debug, Clone, Copy)]
struct Cursor<'a> {
    bytes: &'a [u8],
    position: usize,
}

fn take_u32(cursor: &mut Cursor<'_>) -> Result<u32, Error> {
    let start = cursor.position;
    ensure(start <= cursor.bytes.len(), BAD, "wout file header is damaged")?;
    ensure(cursor.bytes.len() - start >= 4, BAD, "wout file ends inside its header")?;
    let b = cursor.bytes;
    cursor.position = start + 4;
    Ok(u32::from_be_bytes([b[start], b[start + 1], b[start + 2], b[start + 3]]))
}

fn take_count(cursor: &mut Cursor<'_>, limit: usize) -> Result<usize, Error> {
    let count = usize::try_from(take_u32(cursor)?)?;
    ensure(count <= limit, BAD, "wout file header has an oversized entry")?;
    ensure(count <= MAX_FILE_BYTES, BAD, "wout file header is damaged")?;
    Ok(count)
}

fn take_offset(cursor: &mut Cursor<'_>, wide: bool) -> Result<usize, Error> {
    let high = if wide { take_u32(cursor)? } else { 0 };
    let low = take_u32(cursor)?;
    ensure(high == 0, BAD, "wout file is too large for this program")?;
    let offset = usize::try_from(low)?;
    ensure(offset <= cursor.bytes.len(), BAD, "wout file is shorter than its header says")?;
    Ok(offset)
}

fn skip_padded(cursor: &mut Cursor<'_>, length: usize) -> Result<(), Error> {
    ensure(length <= MAX_FILE_BYTES, BAD, "wout file header is damaged")?;
    let padded = length.div_ceil(4) * 4;
    ensure(
        cursor.bytes.len() - cursor.position >= padded,
        BAD,
        "wout file ends inside its header",
    )?;
    cursor.position += padded;
    Ok(())
}

fn take_name<'a>(cursor: &mut Cursor<'a>) -> Result<&'a [u8], Error> {
    ensure(cursor.position <= cursor.bytes.len(), BAD, "wout file header is damaged")?;
    let length = take_count(cursor, MAX_FILE_BYTES)?;
    let start = cursor.position;
    skip_padded(cursor, length)?;
    ensure(length > 0, BAD, "wout file has an entry without a name")?;
    let bytes: &'a [u8] = cursor.bytes;
    Ok(&bytes[start..start + length])
}

fn take_list_header(cursor: &mut Cursor<'_>, tag: u32) -> Result<usize, Error> {
    let found = take_u32(cursor)?;
    let count = take_count(cursor, MAX_NC_ENTRIES)?;
    ensure(
        found == tag || (found == 0 && count == 0),
        BAD,
        "wout file header has an unexpected section",
    )?;
    ensure(count <= MAX_NC_ENTRIES, BAD, "wout file header is damaged")?;
    Ok(count)
}

/// Bytes per value of the number types this program reads; zero for others.
fn value_width(kind: u32) -> usize {
    if kind == TYPE_DOUBLE {
        8
    } else if kind == TYPE_INT || kind == TYPE_FLOAT {
        4
    } else {
        0
    }
}

fn skip_attributes(cursor: &mut Cursor<'_>) -> Result<(), Error> {
    let count = take_list_header(cursor, TAG_ATTRIBUTE)?;
    for _ in 0..count {
        let name = take_name(cursor)?;
        let kind = take_u32(cursor)?;
        let values = take_count(cursor, MAX_FILE_BYTES)?;
        let width =
            value_width(kind) + usize::from(kind == 1 || kind == 2) + 2 * usize::from(kind == 3);
        ensure(width > 0, BAD, "wout file has an attribute of unknown type")?;
        ensure(!name.is_empty(), BAD, "wout file header is damaged")?;
        skip_padded(cursor, values * width)?;
    }
    Ok(())
}

/// Read the file header.
///
/// # Errors
/// Returns an error if the bytes are not a `NetCDF` classic file or the
/// header is damaged or larger than the limits of this program.
pub fn parse_header(bytes: &[u8]) -> Result<Header, Error> {
    ensure(bytes.len() >= 8, BAD, "wout file is too short")?;
    ensure(bytes.len() <= MAX_FILE_BYTES, Code::Limit, "wout file is too large")?;
    ensure(
        bytes[1..4] != *b"HDF",
        BAD,
        "wout file is NetCDF-4 (HDF5); only the NetCDF classic format is supported",
    )?;
    ensure(bytes[0..3] == *b"CDF", BAD, "wout file is not a NetCDF file")?;
    ensure(bytes[3] == 1 || bytes[3] == 2, BAD, "wout file has an unsupported NetCDF version")?;
    let mut cursor = Cursor { bytes, position: 4 };
    let _records = take_u32(&mut cursor)?;
    let n_dims = take_list_header(&mut cursor, TAG_DIMENSION)?;
    ensure(n_dims <= MAX_NC_DIMS, Code::Limit, "wout file has too many dimensions")?;
    let mut dim_len = [0_usize; MAX_NC_DIMS];
    for length in dim_len.iter_mut().take(n_dims) {
        let name = take_name(&mut cursor)?;
        ensure(!name.is_empty(), BAD, "wout file header is damaged")?;
        *length = take_count(&mut cursor, MAX_FILE_BYTES)?;
    }
    skip_attributes(&mut cursor)?;
    let n_variables = take_list_header(&mut cursor, TAG_VARIABLE)?;
    Ok(Header {
        wide_offsets: bytes[3] == 2,
        n_dims,
        dim_len,
        first_variable: cursor.position,
        n_variables,
    })
}

fn take_shape(cursor: &mut Cursor<'_>, header: &Header) -> Result<Variable, Error> {
    let rank = take_count(cursor, MAX_NC_DIMS)?;
    let mut shape = [1_usize; MAX_RANK];
    let mut count = 1_usize;
    let mut fixed = true;
    for axis in 0..rank {
        let id = take_count(cursor, MAX_NC_DIMS)?;
        ensure(id < header.n_dims, BAD, "wout file refers to an unknown dimension")?;
        let length = header.dim_len[id];
        fixed = fixed && length > 0;
        count = count.saturating_mul(length);
        if let Some(entry) = shape.get_mut(axis) {
            *entry = length;
        }
    }
    ensure(count <= MAX_FILE_BYTES, BAD, "wout file header is damaged")?;
    Ok(Variable {
        kind: 0,
        begin: 0,
        // Record variables and variables of too high a rank are marked unusable.
        count: if fixed && rank <= MAX_RANK { count } else { 0 },
        rank,
        shape,
    })
}

/// Locate a variable by name. `None` means the file has no such variable.
///
/// # Errors
/// Returns an error if the variable list is damaged, or the variable is a
/// record variable or has more dimensions than this program supports.
pub fn find_variable(bytes: &[u8], header: &Header, name: &str) -> Result<Option<Variable>, Error> {
    ensure(!name.is_empty(), BAD, "a variable name is required")?;
    ensure(header.first_variable <= bytes.len(), BAD, "wout file header is damaged")?;
    let mut cursor = Cursor { bytes, position: header.first_variable };
    let mut found = None;
    for _ in 0..header.n_variables {
        let candidate = take_name(&mut cursor)?;
        let mut variable = take_shape(&mut cursor, header)?;
        skip_attributes(&mut cursor)?;
        variable.kind = take_u32(&mut cursor)?;
        let _size = take_u32(&mut cursor)?;
        variable.begin = take_offset(&mut cursor, header.wide_offsets)?;
        if candidate == name.as_bytes() && found.is_none() {
            ensure(variable.count > 0, BAD, "a wout variable has an unsupported shape")?;
            found = Some(variable);
        }
    }
    Ok(found)
}

fn value_at(bytes: &[u8], kind: u32, at: usize) -> f64 {
    let b = bytes;
    match kind {
        TYPE_DOUBLE => f64::from_be_bytes([
            b[at],
            b[at + 1],
            b[at + 2],
            b[at + 3],
            b[at + 4],
            b[at + 5],
            b[at + 6],
            b[at + 7],
        ]),
        TYPE_FLOAT => f64::from(f32::from_be_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])),
        _ => f64::from(i32::from_be_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])),
    }
}

/// Read all values of a variable as real numbers into `out`.
///
/// # Errors
/// Returns an error if `out` has the wrong length, the number type is not
/// integer, float or double, or the file is shorter than the variable.
pub fn read_values(bytes: &[u8], variable: &Variable, out: &mut [f64]) -> Result<(), Error> {
    let width = value_width(variable.kind);
    ensure(width > 0, BAD, "a wout variable has an unsupported number type")?;
    ensure(out.len() == variable.count, BAD, "a wout variable does not have the expected size")?;
    ensure(
        variable.begin <= bytes.len() && (bytes.len() - variable.begin) / width >= out.len(),
        BAD,
        "wout file ends inside a variable",
    )?;
    for (index, value) in out.iter_mut().enumerate() {
        *value = value_at(bytes, variable.kind, variable.begin + index * width);
    }
    Ok(())
}
