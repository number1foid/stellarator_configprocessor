//! Allocation guard: the instrument that enforces "no dynamic memory
//! allocation after initialisation".
//!
//! Installed as the global allocator, it forwards every request to the
//! system allocator and counts it. After [`lock`] has been called it refuses
//! every further request, which makes the standard library stop the program.
//! The counters let a test prove that a stretch of code neither requests nor
//! returns heap memory.
//!
//! This is the only code in the repository that contains `unsafe`: the
//! language requires it to implement an allocator. The raw pointers are
//! passed through and never dereferenced.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// Set once start-up is complete; never cleared.
static LOCKED: AtomicBool = AtomicBool::new(false);
/// Number of memory requests seen so far.
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
/// Number of memory blocks returned so far.
static RELEASES: AtomicUsize = AtomicUsize::new(0);

/// The guard. Install it with `#[global_allocator]`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Guard;

// SAFETY: every request is forwarded unchanged to the system allocator, or
// answered with a null pointer, which is the documented way to report that
// no memory is available.
unsafe impl GlobalAlloc for Guard {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let _count = ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        if LOCKED.load(Ordering::Relaxed) {
            return core::ptr::null_mut();
        }
        // SAFETY: the caller upholds the contract of `GlobalAlloc::alloc`.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        let _count = RELEASES.fetch_add(1, Ordering::Relaxed);
        // SAFETY: `block` was obtained from `System.alloc` with `layout`.
        unsafe { System.dealloc(block, layout) }
    }
}

/// Refuse every memory request from now on. Call it when start-up is done.
pub fn lock() {
    LOCKED.store(true, Ordering::Relaxed);
}

/// Whether [`lock`] has been called.
#[must_use]
pub fn is_locked() -> bool {
    LOCKED.load(Ordering::Relaxed)
}

/// Number of memory requests made since the program started.
#[must_use]
pub fn allocations() -> usize {
    ALLOCATIONS.load(Ordering::Relaxed)
}

/// Number of memory blocks returned since the program started.
#[must_use]
pub fn releases() -> usize {
    RELEASES.load(Ordering::Relaxed)
}
