//! Rule 3: no dynamic memory allocation after initialisation.
//!
//! The allocation guard counts every request for and every return of heap
//! memory. This test prepares a complete run, then checks that computing and
//! writing all results changes neither counter.
//!
//! This file holds a single test on purpose: a second test running in
//! another thread would allocate and disturb the counters.

#[path = "support/synthetic.rs"]
mod synthetic;

use std::fs;
use stellarator_configprocessor::cli::{execute, prepare, EXIT_OK};
use stellarator_configprocessor::error::{Code, Error};
use synthetic::{circular_coils_text, torus_wout};

/// Counts memory requests; it is never locked here.
#[global_allocator]
static GUARD: alloc_guard::Guard = alloc_guard::Guard;

#[test]
fn computing_and_writing_allocates_nothing() -> Result<(), Error> {
    let mut folder = std::env::temp_dir();
    folder.push(format!("stellarator_configprocessor_test_{}", std::process::id()));
    fs::create_dir_all(&folder)?;
    let path = |name: &str| folder.join(name).to_string_lossy().into_owned();
    // A circular torus inside 20 circular coils; 0.4 T is the field of the
    // ideal toroidal solenoid on the axis.
    fs::write(path("wout_torus.nc"), torus_wout(10.0, 2.0, 5, 0.4, 11)?)?;
    fs::write(path("coils.torus"), circular_coils_text(20, 10.0, 4.0, 1.0e6, 120, 5)?)?;
    let arguments = [
        path("wout_torus.nc"),
        path("coils.torus"),
        "--wp-radial=0.5".to_owned(),
        "--wp-toroidal=0.5".to_owned(),
        "--epseff=0.01".to_owned(),
        "--neutron-peakfactor=1.6".to_owned(),
        "--points=64".to_owned(),
        "--output".to_owned(),
        path("stella_conf.json"),
        "--report".to_owned(),
        path("report.json"),
    ];
    let mut out = std::io::sink();
    let mut err = std::io::sink();
    let Some(mut session) = prepare(arguments.into_iter(), &mut out)? else {
        return Err(Error::new(Code::Usage, "the command line did not ask for a run"));
    };

    let before = (alloc_guard::allocations(), alloc_guard::releases());
    let status = execute(&mut session, &mut out, &mut err);
    let after = (alloc_guard::allocations(), alloc_guard::releases());

    assert_eq!(status, EXIT_OK, "the run failed");
    assert!(before.0 > 0, "the guard is not counting");
    assert_eq!(after.0 - before.0, 0, "memory was requested during the calculation");
    assert_eq!(after.1 - before.1, 0, "memory was returned during the calculation");
    assert!(!alloc_guard::is_locked(), "the guard must stay unlocked in this test");
    drop(session);
    let config = fs::read_to_string(path("stella_conf.json"))?;
    assert!(config.contains("\"coilspermodule\": 4"), "unexpected configuration: {config}");
    assert!(fs::read_to_string(path("report.json"))?.contains("\"passed\": true"));
    fs::remove_dir_all(&folder)?;
    Ok(())
}
