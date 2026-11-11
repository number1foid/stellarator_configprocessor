//! Checks the product code against "The Power of Ten" rules, as far as they
//! can be decided from the source text. The compiler, Clippy and the
//! allocation test cover the rest; see `docs/power_of_ten.md`.

use std::fs;
use std::path::PathBuf;

/// Longest function allowed by rule 4, in source lines.
const MAX_FUNCTION_LINES: usize = 60;
/// Smallest average number of assertions per function allowed by rule 5.
const MIN_ASSERTION_DENSITY: f64 = 2.0;
/// The only macros the product code may invoke.
const ALLOWED_MACROS: [&str; 3] = ["write", "writeln", "vec"];

/// Constructs that must not appear anywhere in the product code.
const BANNED: [(&str, &str); 34] = [
    ("while ", "rule 2: loops need a fixed upper bound"),
    ("loop ", "rule 2: loops need a fixed upper bound"),
    ("loop{", "rule 2: loops need a fixed upper bound"),
    ("Box<", "rule 3 and 9: no heap indirection"),
    ("Rc<", "rule 3 and 9: no heap indirection"),
    ("Arc<", "rule 3 and 9: no heap indirection"),
    ("RefCell", "rule 9: no hidden indirection"),
    ("unwrap()", "rule 5: a failed check must return an error"),
    (".expect(", "rule 5: a failed check must return an error"),
    ("static mut", "rule 6: no shared mutable data"),
    ("pub static", "rule 6: data at the smallest scope"),
    ("let _ =", "rule 7: a return value may only be discarded by name"),
    ("macro_rules", "rule 8: no macro definitions"),
    ("#[cfg", "rule 8: no conditional compilation"),
    ("#![cfg", "rule 8: no conditional compilation"),
    ("dyn ", "rule 9: no indirect calls"),
    ("fn(", "rule 9: no function pointers"),
    ("impl Fn", "rule 9: no function-valued arguments"),
    ("FnMut", "rule 9: no function-valued arguments"),
    ("FnOnce", "rule 9: no function-valued arguments"),
    ("&mut &", "rule 9: one level of indirection"),
    ("& &", "rule 9: one level of indirection"),
    ("Vec<Vec<", "rule 9: one level of indirection"),
    ("&Vec<", "rule 9: one level of indirection"),
    ("&String", "rule 9: one level of indirection"),
    ("*const ", "rule 9: no raw pointers"),
    (" as ", "rule 10: no lossy casts"),
    ("#[allow", "rule 10: warnings are fixed, not silenced"),
    ("#![allow", "rule 10: warnings are fixed, not silenced"),
    ("#[expect", "rule 10: warnings are fixed, not silenced"),
    ("goto", "rule 1: simple control flow"),
    ("catch_unwind", "rule 1: no non-local jumps"),
    ("std::thread", "rule 3: threads allocate"),
    ("std::process::exit", "rule 1: one way out of the program"),
];

/// One file of the product code.
struct Source {
    name: String,
    /// The text without comments and without the contents of literals.
    code: String,
}

/// One function of the product code.
struct Function {
    name: String,
    file: String,
    lines: usize,
    body: String,
}

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Remove comments and blank out the contents of string and character
/// literals, keeping every line break.
fn strip(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::new();
    let mut index = 0;
    while index < chars.len() {
        let (here, next) = (chars[index], chars.get(index + 1).copied().unwrap_or(' '));
        if here == '/' && next == '/' {
            while index < chars.len() && chars[index] != '\n' {
                index += 1;
            }
        } else if here == '"' {
            out.push('"');
            index += 1;
            while index < chars.len() && chars[index] != '"' {
                if chars[index] == '\n' {
                    out.push('\n');
                }
                index += if chars[index] == '\\' { 2 } else { 1 };
            }
            out.push('"');
            index += 1;
        } else if here == '\'' && next == '\\' {
            out.push_str("' '");
            index += 2;
            while index < chars.len() && chars[index] != '\'' {
                index += 1;
            }
            index += 1;
        } else if here == '\'' && chars.get(index + 2) == Some(&'\'') {
            out.push_str("' '");
            index += 3;
        } else {
            out.push(here);
            index += 1;
        }
    }
    out
}

/// The product code: the library, the program and the allocation guard.
fn sources() -> Vec<Source> {
    let mut paths = vec![root().join("alloc_guard").join("src").join("lib.rs")];
    for entry in fs::read_dir(root().join("src")).into_iter().flatten().flatten() {
        if entry.path().extension().is_some_and(|extension| extension == "rs") {
            paths.push(entry.path());
        }
    }
    paths.sort();
    let mut found = Vec::new();
    for path in paths {
        let text = fs::read_to_string(&path).unwrap_or_default();
        assert!(!text.is_empty(), "cannot read {}", path.display());
        let guard = path.components().any(|part| part.as_os_str() == "alloc_guard");
        let file = path.file_name().map(|name| name.to_string_lossy().into_owned());
        let name =
            format!("{}{}", if guard { "alloc_guard/" } else { "" }, file.unwrap_or_default());
        found.push(Source { name, code: strip(&text) });
    }
    assert!(found.len() >= 15, "expected the product sources, found {} files", found.len());
    found
}

/// The name defined by a line that starts a function, if it does.
fn function_name(line: &str) -> Option<&str> {
    let mut rest = line.trim_start();
    for prefix in ["pub(crate) ", "pub ", "const ", "unsafe "] {
        rest = rest.strip_prefix(prefix).unwrap_or(rest);
    }
    let rest = rest.strip_prefix("fn ")?;
    let end = rest.find(|character: char| !(character.is_alphanumeric() || character == '_'))?;
    Some(&rest[..end])
}

/// Every function of the product code. The code is formatted by rustfmt, so
/// a function ends at the first closing brace with the indentation of its
/// first line.
fn functions(sources: &[Source]) -> Vec<Function> {
    let mut found = Vec::new();
    for source in sources {
        let lines: Vec<&str> = source.code.lines().collect();
        for (start, line) in lines.iter().enumerate() {
            let Some(name) = function_name(line) else { continue };
            let indent = line.len() - line.trim_start().len();
            let closing = format!("{}}}", " ".repeat(indent));
            let length = lines[start..].iter().position(|other| *other == closing);
            let length = length.map_or(1, |offset| offset + 1);
            found.push(Function {
                name: name.to_owned(),
                file: source.name.clone(),
                lines: length,
                body: lines[start..start + length].join("\n"),
            });
        }
    }
    assert!(found.len() >= 100, "expected the product functions, found {}", found.len());
    found
}

/// Number of calls of `name` in `code`: the name followed by an opening
/// parenthesis and not preceded by a letter, digit or underscore.
fn count_calls(code: &str, name: &str) -> usize {
    let pattern = format!("{name}(");
    let mut count = 0;
    for (position, _) in code.match_indices(&pattern) {
        let before = code[..position].chars().next_back().unwrap_or(' ');
        // A call written with a dot is a method of a standard type: the
        // product code defines no methods of its own (checked below).
        if !(before.is_alphanumeric() || before == '_' || before == '.') {
            count += 1;
        }
    }
    count
}

#[test]
fn rule_4_functions_fit_on_one_page() {
    for function in functions(&sources()) {
        assert!(
            function.lines <= MAX_FUNCTION_LINES,
            "{} in {} has {} lines; the limit is {MAX_FUNCTION_LINES}",
            function.name,
            function.file,
            function.lines
        );
    }
}

#[test]
fn rule_5_assertion_density() {
    let sources = sources();
    let count = functions(&sources).len();
    let mut assertions = 0;
    for source in &sources {
        assertions += count_calls(&source.code, "ensure") + count_calls(&source.code, "ensure_at");
        // The definitions of the two assertion functions are not assertions.
        if source.name == "error.rs" {
            assertions -= 2;
        }
        for banned in ["assert!", "assert_eq!", "debug_assert", "panic!", "unreachable!", "todo!"] {
            assert!(!source.code.contains(banned), "{}: `{banned}` is not allowed", source.name);
        }
    }
    let density = f64::from(u32::try_from(assertions).unwrap_or(0))
        / f64::from(u32::try_from(count).unwrap_or(1));
    assert!(
        density >= MIN_ASSERTION_DENSITY,
        "{assertions} assertions in {count} functions: {density:.3} per function; the minimum \
         is {MIN_ASSERTION_DENSITY}"
    );
}

#[test]
fn rules_1_to_10_banned_constructs() {
    for source in sources() {
        let guard = source.name.starts_with("alloc_guard/");
        for (token, rule) in BANNED {
            assert!(!source.code.contains(token), "{}: `{token}` ({rule})", source.name);
        }
        if !guard {
            for token in ["unsafe", "*mut "] {
                assert!(
                    !source.code.contains(token),
                    "{}: `{token}` is only allowed in the allocation guard",
                    source.name
                );
            }
            for line in source.code.lines() {
                let item = line.trim_start().starts_with("static ");
                assert!(
                    !item || source.name == "main.rs",
                    "{}: a static item is only allowed for the allocator (rule 6)",
                    source.name
                );
            }
        }
        let chars: Vec<char> = source.code.chars().collect();
        for (index, character) in chars.iter().enumerate() {
            let before = if index > 0 { chars[index - 1] } else { ' ' };
            let after = chars.get(index + 1).copied().unwrap_or(' ');
            assert!(
                *character != '|' || before == '|' || after == '|',
                "{}: closure or pattern alternative (rule 9: no function-valued arguments)",
                source.name
            );
            assert!(
                !(*character == '&' && before == '&' && after != ' ' && after != '\n'),
                "{}: reference to a reference (rule 9: one level of indirection)",
                source.name
            );
        }
    }
}

#[test]
fn rule_8_only_standard_macros() {
    for source in sources() {
        let chars: Vec<char> = source.code.chars().collect();
        for (index, character) in chars.iter().enumerate() {
            let opens = matches!(chars.get(index + 1), Some('(' | '[' | '{'));
            if *character != '!' || !opens || index == 0 {
                continue;
            }
            let mut start = index;
            while start > 0 && (chars[start - 1].is_alphanumeric() || chars[start - 1] == '_') {
                start -= 1;
            }
            let name: String = chars[start..index].iter().collect();
            assert!(
                name.is_empty() || ALLOWED_MACROS.contains(&name.as_str()),
                "{}: macro `{name}!` is not allowed (rule 8)",
                source.name
            );
        }
    }
}

#[test]
fn rule_2_loops_over_input_are_capped() {
    for source in sources() {
        for line in source.code.lines() {
            let over_input = line.contains(".lines()") || line.contains(".chars()");
            if line.trim_start().starts_with("for ") && (over_input || line.contains("arguments")) {
                assert!(line.contains(".take("), "{}: uncapped loop: {}", source.name, line.trim());
            }
        }
    }
}

#[test]
fn rule_1_no_recursion() {
    let all = functions(&sources());
    for function in &all {
        let signature = function.body.split_once('{').map_or("", |parts| parts.0);
        let method = signature.contains("self,") || signature.contains("self)");
        let required = ["fmt", "alloc", "dealloc"].contains(&function.name.as_str());
        assert!(
            !method || required,
            "{} is a method; only free functions are allowed",
            function.name
        );
    }
    let mut names: Vec<&str> = all.iter().map(|function| function.name.as_str()).collect();
    names.sort_unstable();
    names.dedup();
    // Functions with the same name are merged, which can only add edges.
    let mut callees: Vec<Vec<usize>> = vec![Vec::new(); names.len()];
    for function in &all {
        let caller = names.binary_search(&function.name.as_str()).unwrap_or(0);
        let body = function.body.split_once('{').map_or("", |parts| parts.1);
        for (callee, name) in names.iter().enumerate() {
            if count_calls(body, name) > 0 && !callees[caller].contains(&callee) {
                callees[caller].push(callee);
            }
        }
    }
    // Repeatedly remove functions that call nothing that is left. In a call
    // graph without cycles every function is removed.
    let mut left = vec![true; names.len()];
    for _ in 0..names.len() {
        for index in 0..names.len() {
            if left[index] && !callees[index].iter().any(|callee| left[*callee]) {
                left[index] = false;
            }
        }
    }
    let cyclic: Vec<&str> =
        names.iter().zip(left.iter()).filter(|pair| *pair.1).map(|pair| *pair.0).collect();
    assert!(cyclic.is_empty(), "functions on a call cycle (rule 1): {cyclic:?}");
}

#[test]
fn rules_1_3_and_10_build_settings() {
    let manifest = fs::read_to_string(root().join("Cargo.toml")).unwrap_or_default();
    assert_eq!(manifest.matches("panic = \"abort\"").count(), 2, "rule 1: panics must abort");
    for setting in [
        "warnings = { level = \"deny\"",
        "pedantic = { level = \"deny\"",
        "all = { level = \"deny\"",
        "unsafe_code = \"forbid\"",
        "unused_results = \"deny\"",
        "as_conversions = \"deny\"",
    ] {
        assert!(manifest.contains(setting), "Cargo.toml lacks `{setting}` (rule 10)");
    }
    let clippy = fs::read_to_string(root().join("clippy.toml")).unwrap_or_default();
    assert!(clippy.contains("too-many-lines-threshold = 60"), "rule 4: Clippy line limit");
    let program = fs::read_to_string(root().join("src").join("main.rs")).unwrap_or_default();
    assert!(program.contains("#[global_allocator]"), "rule 3: the guard must be installed");
    assert!(program.contains("alloc_guard::lock()"), "rule 3: the allocator must be locked");
}
