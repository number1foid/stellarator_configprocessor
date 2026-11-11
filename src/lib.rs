//! Pre-processor that writes the PROCESS stellarator configuration file.
//!
//! The program reads a VMEC `wout` equilibrium and a MAKEGRID coil file and
//! computes the values that the PROCESS systems code reads from
//! `stella_conf.json` when `istell = 6`.
//!
//! The code follows the ten rules of "The Power of Ten - Rules for
//! Developing Safety Critical Code" (G. J. Holzmann, NASA/JPL). How each rule
//! is met, and how that is checked, is described in `docs/power_of_ten.md`.

pub mod cli;
pub mod coilgeom;
pub mod coils;
pub mod consistency;
pub mod curve;
pub mod energy;
pub mod equilibrium;
pub mod error;
pub mod field;
pub mod forces;
pub mod limits;
pub mod makegrid;
pub mod netcdf;
pub mod num;
pub mod pipeline;
pub mod plasma;
pub mod report;
pub mod schema;
pub mod winding_pack;
