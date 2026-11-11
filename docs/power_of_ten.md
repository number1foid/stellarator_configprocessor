# Compliance with "The Power of Ten"

This document states how the Rust code in this repository meets each of the
ten rules of G. J. Holzmann, *The Power of Ten - Rules for Developing Safety
Critical Code* (NASA/JPL Laboratory for Reliable Software), and how each rule
is checked.

The rules were written for C. Where a rule names a C construct, the nearest
Rust construct is treated the same way; those readings are marked
**Interpretation** below so that they can be judged.

## Scope

| Code | Covered | Remark |
|---|---|---|
| `src/` (library and program) | yes | All ten rules, all checks below |
| `alloc_guard/` (allocation guard) | yes, with one stated exception | Contains the only `unsafe` code; see rule 9 |
| `tests/` | no | Test code is formatted and linted at the same level, but it allocates freely, uses `assert!` and closures |
| `scripts/` (Python helpers) | no | They fetch public test data and drive PROCESS, which is itself Python; they are not part of the program |
| The Rust standard library | no | As with the C library in the original paper, it is used but not written to these rules |

No third-party crate is used. None is written to these rules, and any
dependency would have fallen outside what can be claimed here. The command
line parser, the JSON writer and the reader for the NetCDF classic format are
therefore part of this code base.

## Figures

Measured by `tests/power_of_ten.rs` on the product code:

| Quantity | Value | Limit |
|---|---|---|
| Functions | 176 | |
| Longest function | 60 lines | 60 |
| Assertions (`ensure`, `ensure_at`) | 355 | |
| Assertions per function | 2.02 | at least 2 |
| `unsafe` blocks outside the allocation guard | 0 | 0 |
| Compiler and Clippy findings | 0 | 0 |
| Heap requests during the calculation | 0 | 0 |

## The rules

### 1. Simple control flow; no goto, setjmp/longjmp or recursion

- The call graph has no cycle. `rule_1_no_recursion` builds the call graph
  of all product functions from the source and fails if any function can
  reach itself. Functions that share a name are merged, which can only add
  edges, so the check errs on the safe side.
- Rust has no `goto`. Its counterpart of `setjmp`/`longjmp` is unwinding
  after a panic. Both build profiles set `panic = "abort"`, so a failed
  run-time check of the language (for example an index out of range) stops
  the program at once and nothing is unwound. `catch_unwind` and
  `std::process::exit` are banned; the program has one way out, the return
  value of `main`.
- The product code contains no `panic!`, `unwrap` or `expect`; Clippy
  rejects them.

### 2. All loops have a fixed upper bound

- Every loop is a `for` loop over a range or a slice. `while` and `loop` are
  banned by `rules_1_to_10_banned_constructs`.
- The length of every range and slice is checked against a constant in
  [`src/limits.rs`](../src/limits.rs) when the data is read or the buffer is
  created. Loops driven by the content of an input file or the command line
  are capped with `.take(LIMIT)`; `rule_2_loops_over_input_are_capped`
  checks that.
- Searches that would naturally be open-ended (bisection) run a fixed number
  of steps (`BISECTION_STEPS`).

**Interpretation and limit.** The rule asks that a tool can prove each bound
statically. No Rust tool does that today. What is checked mechanically is
the form of the loops (bounded constructs only); that every length traces
back to a constant was established by review, and is enforced at run time by
the checks at the points where data enters.

### 3. No dynamic memory allocation after initialisation

- The program runs in two phases. `cli::prepare` reads the command line and
  both input files, allocates every buffer (`pipeline::new_workspace`) and
  opens the output files. `cli::execute` computes and writes the results.
- Between the two, `main` locks the allocator (`alloc_guard::lock`). From
  then on any request for memory stops the program. Every run of the
  program therefore enforces the rule.
- `tests/no_allocation.rs` runs a complete case and checks that the number
  of memory requests and of memory returns does not change between the
  start and the end of `execute`, including writing the configuration file
  and the JSON report.
- Errors carry no heap data, so they can be created and reported after the
  lock.

**Interpretation.** "Initialisation" includes reading the two input files,
because the sizes of the buffers follow from them. All sizes are bounded by
the constants in `limits.rs`. Memory is returned to the system only when the
program ends.

### 4. No function longer than 60 lines

- `rule_4_functions_fit_on_one_page` measures every function from its first
  line to its closing brace, signature, comments and blank lines included,
  and fails above 60 lines.
- Clippy's `too_many_lines` is set to the same limit in `clippy.toml`.
- The code is formatted by `rustfmt` with one statement per line; CI fails
  if it is not.

### 5. At least two assertions per function on average

- The assertion of this code base is `ensure(condition, code, message)?`
  (and `ensure_at`, which adds a line number). It is a side-effect-free
  Boolean test. When it fails, the function returns an error to its caller,
  which is the recovery action the rule asks for; the error travels up to
  `main`, which reports it and returns a non-zero exit code.
- `rule_5_assertion_density` counts the assertions and the functions and
  fails below an average of 2. `assert!`, `debug_assert!` and `panic!` are
  not used: they abort instead of recovering.
- The assertions check parameters, sizes of buffers against each other,
  indices, and that results are finite and in range. None is of the form
  "true is true"; each one can fail for some input.

**Limit.** Small arithmetic helpers (vector sum, scalar product and the
like) have no assertion; the average is carried by the functions that call
them. No tool proves that every assertion can fail; that was established by
review.

### 6. Data declared at the smallest possible scope

- There is no global mutable data: `static mut` and `pub static` are
  banned. The only `static` item is the allocator in `main.rs`, which the
  language requires to be global, and the three private statics of the allocation guard.
- All working memory is owned by one `Workspace` value created in
  `prepare` and passed down explicitly. Functions receive only the parts
  they use.
- Variables are immutable unless declared `mut`; the compiler rejects any
  variable that is declared `mut` without need or never used.

### 7. Return values checked; parameters validated

- The compiler lints `unused_must_use` and `unused_results` are errors, so
  no return value can be dropped silently. A value may be discarded only by
  binding it to a name that starts with an underscore, which states the
  intent; the anonymous `let _ =` is banned.
- Every fallible function returns `Result`, and every call site either
  handles the error or passes it on with `?`.
- Parameters are validated inside each function by the assertions of
  rule 5.

### 8. Limited use of the preprocessor

- Rust has no textual preprocessor. Its counterpart is the macro system and
  conditional compilation.
- The product code defines no macro (`macro_rules` is banned) and contains
  no conditional compilation (`#[cfg` is banned), so exactly one version of
  the code exists.
- It invokes three macros of the standard library: `write!`, `writeln!` and
  `vec!` (the last only during initialisation).
  `rule_8_only_standard_macros` fails on any other macro.

**Interpretation.** Invoking `write!` is treated like calling `printf` in
the paper's own example. `#[derive(...)]` on data types is used; it is
compiler-generated code for copying and printing values.

### 9. Restricted use of pointers; no function pointers

- `unsafe_code = "forbid"` for the library and the program: no raw pointer
  can be created or dereferenced there. References are checked by the
  compiler for validity and aliasing.
- No function pointers, no trait objects (`dyn`), no closures and no
  function-valued arguments. Every call is a direct call to a named
  function, so the call graph used for rule 1 is complete.
- No heap indirection types (`Box`, `Rc`, `Arc`, `RefCell`), no reference
  to a reference, no nested containers (`Vec<Vec<..>>`). Arrays are flat and
  indexed arithmetically.

**Interpretation.** "One level of dereferencing" is read as one level of
indirection in every type. A function may receive a reference to a record
whose fields own flat arrays, the Rust form of `p->array[i]`. Output
functions are generic over the writer (`W: Write`), which the compiler
resolves to direct calls.

**Exception.** `alloc_guard` implements the allocator interface, which the
language only allows in `unsafe` code. It passes raw pointers through to the
system allocator and never dereferences them.

### 10. All warnings enabled; static analysis with zero findings

- `Cargo.toml` turns every compiler warning into an error and enables
  further lints that are off by default (missing documentation, unused
  results, trivial casts, unreachable public items).
- Clippy runs with the groups `all` and `pedantic` as errors, plus
  restriction lints against `unwrap`, `expect`, `panic`, `as` casts,
  printing macros and `exit`.
- No finding is silenced: `#[allow]` and `#[expect]` are banned in the
  product code. Where a lint objected, the code was rewritten.
- Three analysers run on every push: the compiler, Clippy, and the rule
  checker in `tests/power_of_ten.rs`. `rustfmt --check` runs as well.

**Limit.** The workflow also defines a daily run, but GitHub starts
scheduled runs only from the default branch. Until this branch is the
default, the analysis runs on every push and not daily.

## How to run the checks

```
cargo fmt --all -- --check
cargo clippy --workspace --all-targets
cargo test --workspace --release
```

`cargo test` runs the rule checker, the allocation test and the functional
tests. To see the rule checker alone: `cargo test --release --test
power_of_ten`.
