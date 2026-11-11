//! Fixed upper bounds on every size in the program.
//!
//! Every loop runs over a range or a slice whose length was checked against
//! one of these constants when the data was read or the buffers were created,
//! so the number of iterations of every loop has a compile-time bound.

/// Largest size of an input file [bytes].
pub const MAX_FILE_BYTES: usize = 1 << 30;
/// Largest number of lines read from a coil file.
pub const MAX_LINES: usize = 4_000_000;
/// Largest number of whitespace-separated fields inspected on one line.
pub const MAX_TOKENS: usize = 64;
/// Largest number of characters of a number in a coil file.
pub const MAX_NUMBER_CHARS: usize = 64;
/// Largest number of coils in the full torus.
pub const MAX_COILS: usize = 512;
/// Largest number of points of one coil in the coil file.
pub const MAX_FILE_POINTS: usize = 8192;
/// Largest number of samples per coil.
pub const MAX_SAMPLES: usize = 1024;
/// Smallest number of samples per coil.
pub const MIN_SAMPLES: usize = 32;
/// Points used for a one-line (circular) coil in the coil file.
pub const CIRCLE_POINTS: usize = 128;
/// Largest number of Fourier modes of the equilibrium.
pub const MAX_MODES: usize = 8192;
/// Largest number of radial grid points of the equilibrium.
pub const MAX_RADIAL: usize = 4096;
/// Largest number of field periods.
pub const MAX_PERIODS: usize = 1000;
/// Largest number of filaments across one side of the winding pack.
pub const MAX_SUB: usize = 12;
/// Largest number of field points across one side of the winding pack.
pub const MAX_CORNERS: usize = MAX_SUB + 1;
/// Largest number of winding pack sizes used for the peak-field fit.
pub const MAX_SCAN: usize = 8;
/// Largest number of surface points per angle.
pub const MAX_SURFACE: usize = 512;
/// Number of points along the magnetic axis for the consistency check.
pub const AXIS_POINTS: usize = 720;
/// Flux surfaces used for the derivative of the plasma-coil distance.
pub const FIT_SURFACES: usize = 4;
/// Largest number of samples of a coil for the self-inductance integral.
pub const MAX_INDUCTANCE_SAMPLES: usize = 16_384;
/// Largest number of `GROUP=FACTOR` settings on the command line.
pub const MAX_GROUPS: usize = 32;
/// Largest number of command-line arguments.
pub const MAX_ARGS: usize = 128;
/// Largest number of dimensions declared in a `NetCDF` file.
pub const MAX_NC_DIMS: usize = 256;
/// Largest number of variables or attributes declared in a `NetCDF` file.
pub const MAX_NC_ENTRIES: usize = 16_384;
/// Largest rank of a variable this program reads.
pub const MAX_RANK: usize = 4;
/// Steps of a bisection over a 64-bit range.
pub const BISECTION_STEPS: usize = 64;
/// Largest number of characters of a path or name written to a report.
pub const MAX_TEXT: usize = 4096;
