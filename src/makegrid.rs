//! Reader for MAKEGRID `coils.*` files.
//!
//! File layout (after an optional Fortran namelist header):
//!
//! ```text
//! periods 5
//! begin filament
//! mirror NIL
//! x y z current
//! ...
//! x y z 0 group_id group_name      <- closes the coil
//! R 0 z current group_id name      <- one-line coil: circle around the z axis
//! end
//! ```
//!
//! The file lists every coil of the full torus. The current of a coil is the
//! value in the fourth column, optionally multiplied by a per-group factor
//! (the `extcur` array of VMEC). Coils without current are left out.

use crate::error::{ensure, ensure_at, Code, Error};
use crate::limits::{
    CIRCLE_POINTS, MAX_COILS, MAX_FILE_BYTES, MAX_FILE_POINTS, MAX_GROUPS, MAX_LINES,
    MAX_NUMBER_CHARS, MAX_PERIODS, MAX_TOKENS,
};
use crate::num::{finite_vector, floor_count, norm, real, sub, Vec3};

const BAD: Code = Code::CoilFile;

/// Factors applied to the currents of coil groups (`GROUP=FACTOR`).
#[derive(Debug, Clone, Copy)]
pub struct GroupFactors {
    /// Number of entries in use.
    pub count: usize,
    /// Group identifiers.
    pub group: [f64; MAX_GROUPS],
    /// Factor for each group.
    pub factor: [f64; MAX_GROUPS],
}

/// Irregularities met while reading the file; reported as warnings.
#[derive(Debug, Clone, Copy, Default)]
pub struct FileNotes {
    /// Filament coils left out because their current is zero.
    pub zero_current: usize,
    /// Coils whose closing line carries a non-zero current.
    pub closing_current: usize,
    /// Coils whose current varies along the filament (first value is used).
    pub varying_current: usize,
    /// Coils that are not closed (closed with a straight segment).
    pub not_closed: usize,
}

/// All current-carrying coils of the full torus, as read from the file.
#[derive(Debug, Clone)]
pub struct CoilSet {
    /// Number of field periods declared in the file.
    pub periods: usize,
    /// Index of the first point of each coil, plus one final entry.
    pub offsets: Vec<usize>,
    /// Distinct points of all coils; each curve closes implicitly `[m]`.
    pub points: Vec<Vec3>,
    /// Current of each coil `[A]` (ampere-turns).
    pub currents: Vec<f64>,
    /// Irregularities met while reading.
    pub notes: FileNotes,
}

/// What one line of the file contains.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Skip,
    End,
    Point,
    Closing,
}

/// One parsed line.
#[derive(Debug, Clone, Copy)]
struct Line {
    kind: Kind,
    point: Vec3,
    current: f64,
    group: f64,
}

/// The coil being read.
#[derive(Debug, Clone, Copy)]
struct Pending {
    begin: usize,
    first: f64,
    low: f64,
    high: f64,
}

/// A number in Fortran or C notation (`1.5D+03` or `1.5E+03`).
fn parse_real(token: &str, line: usize) -> Result<f64, Error> {
    ensure_at(!token.is_empty(), BAD, "a number is missing", line)?;
    ensure_at(token.len() <= MAX_NUMBER_CHARS, BAD, "a number is too long", line)?;
    let mut buffer = [b' '; MAX_NUMBER_CHARS];
    for (slot, byte) in buffer.iter_mut().zip(token.bytes()) {
        *slot = if byte == b'D' || byte == b'd' { b'E' } else { byte };
    }
    let Ok(text) = core::str::from_utf8(&buffer[..token.len()]) else {
        return Err(Error { code: BAD, message: "cannot parse numbers", line });
    };
    let Ok(value) = text.parse::<f64>() else {
        return Err(Error { code: BAD, message: "cannot parse numbers", line });
    };
    ensure_at(value.is_finite(), BAD, "non-finite number", line)?;
    Ok(value)
}

fn find_periods(text: &str) -> Result<(usize, usize), Error> {
    ensure(text.len() <= MAX_FILE_BYTES, Code::Limit, "the coil file is too large")?;
    let mut found = None;
    for (index, line) in text.lines().enumerate().take(MAX_LINES) {
        let mut tokens = line.split_whitespace();
        let first = tokens.next().unwrap_or("");
        if found.is_none() && first.eq_ignore_ascii_case("periods") {
            let value = parse_real(tokens.next().unwrap_or(""), index + 1)?;
            ensure_at(value >= 1.0 && value.fract() == 0.0, BAD, "bad 'periods' line", index + 1)?;
            found = Some((floor_count(value, MAX_PERIODS)?, index + 1));
        }
    }
    ensure(found.is_some(), BAD, "no 'periods' line found; not a MAKEGRID coils file")?;
    found.ok_or(Error::new(BAD, "no 'periods' line found; not a MAKEGRID coils file"))
}

fn parse_line(text: &str, number: usize) -> Result<Line, Error> {
    let mut line = Line { kind: Kind::Skip, point: [0.0; 3], current: 0.0, group: 0.0 };
    let mut tokens = text.split_whitespace().take(MAX_TOKENS);
    let first = tokens.next().unwrap_or("");
    if first.eq_ignore_ascii_case("end") {
        line.kind = Kind::End;
    }
    let keyword = first.is_empty()
        || first.eq_ignore_ascii_case("begin")
        || first.eq_ignore_ascii_case("mirror")
        || line.kind == Kind::End;
    if keyword {
        return Ok(line);
    }
    line.point[0] = parse_real(first, number)?;
    let (y, z, current) = (tokens.next(), tokens.next(), tokens.next());
    ensure_at(current.is_some(), BAD, "expected 'x y z current'", number)?;
    line.point[1] = parse_real(y.unwrap_or(""), number)?;
    line.point[2] = parse_real(z.unwrap_or(""), number)?;
    line.current = parse_real(current.unwrap_or(""), number)?;
    line.kind = Kind::Point;
    if let Some(group) = tokens.next() {
        let Ok(value) = group.parse::<f64>() else {
            return Err(Error { code: BAD, message: "bad group id", line: number });
        };
        ensure_at(value.is_finite(), BAD, "bad group id", number)?;
        line.group = value.trunc();
        line.kind = Kind::Closing;
    }
    Ok(line)
}

fn group_factor(factors: &GroupFactors, group: f64) -> Result<f64, Error> {
    ensure(factors.count <= MAX_GROUPS, Code::Limit, "too many group factors")?;
    ensure(group.is_finite(), BAD, "bad group id")?;
    let mut factor = 1.0;
    for index in 0..factors.count {
        if (factors.group[index] - group).abs() < 0.5 {
            factor = factors.factor[index];
        }
    }
    Ok(factor)
}

/// One-line coil: a circle around the z axis through the given point.
fn push_circle(set: &mut CoilSet, line: &Line, current: f64, number: usize) -> Result<(), Error> {
    let radius = line.point[0].hypot(line.point[1]);
    ensure_at(radius > 0.0, BAD, "a circular coil needs a radius", number)?;
    ensure_at(set.currents.len() < MAX_COILS, Code::Limit, "too many coils", number)?;
    for index in 0..CIRCLE_POINTS {
        let angle = core::f64::consts::TAU * real(index)? / real(CIRCLE_POINTS)?;
        set.points.push([radius * angle.cos(), radius * angle.sin(), line.point[2]]);
    }
    set.offsets.push(set.points.len());
    set.currents.push(current);
    Ok(())
}

/// Remove repeated points from the coil that starts at `begin`.
///
/// Consecutive duplicates and a final point that repeats the first are
/// dropped. Returns whether the coil was open in the file.
fn clean_points(set: &mut CoilSet, begin: usize, number: usize) -> Result<bool, Error> {
    ensure_at(begin < set.points.len(), BAD, "a coil has no points", number)?;
    let tail = &mut set.points[begin..];
    let (mut low, mut high) = (tail[0], tail[0]);
    for point in tail.iter() {
        for axis in 0..3 {
            low[axis] = low[axis].min(point[axis]);
            high[axis] = high[axis].max(point[axis]);
        }
    }
    let extent = (high[0] - low[0]).max(high[1] - low[1]).max(high[2] - low[2]);
    ensure_at(extent > 0.0, BAD, "curve has no extent", number)?;
    let open = norm(sub(tail[tail.len() - 1], tail[0])) > 1.0e-6 * extent;
    let tolerance = 1.0e-9 * extent;
    let mut kept = 1;
    let mut previous = tail[0];
    for index in 1..tail.len() {
        let point = tail[index];
        if norm(sub(point, previous)) > tolerance {
            tail[kept] = point;
            kept += 1;
        }
        previous = point;
    }
    if kept > 1 && norm(sub(tail[kept - 1], tail[0])) <= tolerance {
        kept -= 1;
    }
    ensure_at(kept >= 3, BAD, "a closed curve needs at least 3 distinct points", number)?;
    set.points.truncate(begin + kept);
    Ok(open)
}

/// Close the coil being read: record irregularities and store it.
fn finish_coil(
    set: &mut CoilSet,
    pending: &Pending,
    closing: &Line,
    factor: f64,
    number: usize,
) -> Result<(), Error> {
    ensure_at(set.currents.len() < MAX_COILS, Code::Limit, "too many coils", number)?;
    ensure_at(factor.is_finite(), BAD, "a group factor is not finite", number)?;
    if closing.current != 0.0 {
        set.notes.closing_current += 1;
    }
    let largest = pending.low.abs().max(pending.high.abs()).max(1.0);
    if pending.high - pending.low > 1.0e-9 * largest {
        set.notes.varying_current += 1;
    }
    if clean_points(set, pending.begin, number)? {
        set.notes.not_closed += 1;
    }
    let current = pending.first * factor;
    if current == 0.0 {
        set.points.truncate(pending.begin);
        set.notes.zero_current += 1;
    } else {
        set.offsets.push(set.points.len());
        set.currents.push(current);
    }
    Ok(())
}

fn add_point(
    set: &mut CoilSet,
    pending: &mut Option<Pending>,
    line: &Line,
    number: usize,
) -> Result<(), Error> {
    ensure_at(finite_vector(line.point), BAD, "non-finite number", number)?;
    let begin = set.points.len();
    let mut state = pending.unwrap_or(Pending {
        begin,
        first: line.current,
        low: line.current,
        high: line.current,
    });
    ensure_at(
        begin - state.begin < MAX_FILE_POINTS,
        Code::Limit,
        "a coil has too many points",
        number,
    )?;
    state.low = state.low.min(line.current);
    state.high = state.high.max(line.current);
    set.points.push(line.point);
    *pending = Some(state);
    Ok(())
}

/// Read a MAKEGRID coil file held in memory.
///
/// # Errors
/// Returns an error if the text is not a MAKEGRID coil file, a line cannot
/// be parsed, a coil has fewer than three distinct points, the file ends
/// inside a coil, no coil carries current, or a size exceeds the limits.
pub fn read_makegrid(text: &str, factors: &GroupFactors) -> Result<CoilSet, Error> {
    let (periods, start) = find_periods(text)?;
    let mut set = CoilSet {
        periods,
        offsets: vec![0],
        points: Vec::new(),
        currents: Vec::new(),
        notes: FileNotes::default(),
    };
    let mut pending: Option<Pending> = None;
    let mut ended = false;
    for (index, text_line) in text.lines().enumerate().take(MAX_LINES).skip(start) {
        let number = index + 1;
        let line = if ended { parse_line("", number)? } else { parse_line(text_line, number)? };
        ended = ended || line.kind == Kind::End;
        if line.kind == Kind::Point {
            add_point(&mut set, &mut pending, &line, number)?;
        }
        if line.kind == Kind::Closing {
            let factor = group_factor(factors, line.group)?;
            if let Some(state) = pending {
                set.points.push(line.point);
                finish_coil(&mut set, &state, &line, factor, number)?;
            } else if line.current * factor != 0.0 {
                push_circle(&mut set, &line, line.current * factor, number)?;
            }
            pending = None;
        }
    }
    ensure(pending.is_none(), BAD, "file ends inside a coil (no closing line with group id)")?;
    ensure(!set.currents.is_empty() || set.notes.zero_current > 0, BAD, "no coils found")?;
    ensure(!set.currents.is_empty(), BAD, "no coil carries current")?;
    Ok(set)
}

/// Multiply all lengths and all currents of the coil set by the factors.
///
/// # Errors
/// Returns an error if a factor is not finite or the length factor is not
/// positive.
pub fn scale_coils(set: &mut CoilSet, length: f64, current: f64) -> Result<(), Error> {
    ensure(
        length.is_finite() && length > 0.0,
        Code::Numeric,
        "the length factor must be positive",
    )?;
    ensure(
        current.is_finite() && current != 0.0,
        Code::Numeric,
        "the current factor must not be zero",
    )?;
    for point in &mut set.points {
        point[0] *= length;
        point[1] *= length;
        point[2] *= length;
    }
    for value in &mut set.currents {
        *value *= current;
    }
    Ok(())
}
