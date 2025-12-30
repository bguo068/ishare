# Error Handling Improvement Plan for ishare

## Executive Summary

This plan outlines improvements to error handling using snafu to enhance structure, reduce memory overhead, and better support context and backtraces.

## Current State Analysis

### Statistics
- **47 files** with Snafu-based error types
- **99 occurrences** of `Option<Backtrace>` fields
- **2 occurrences** of non-optional `Backtrace` fields
- Current snafu version: **0.8.9** (latest stable)

### Current Patterns Observed

1. **Backtrace Inconsistency**: Most errors use `backtrace: Option<Backtrace>`, which adds memory overhead even when backtraces aren't captured
2. **Mixed Error Structures**: Three different patterns:
   - Simple errors without context (e.g., `io.rs`)
   - Errors with rich context and custom display messages (e.g., `share/ibd.rs`)
   - Errors with some context fields (e.g., `gmap.rs`)
3. **Inconsistent Display Messages**: Some errors use `#[snafu(display())]` for clear messages, others rely on default formatting
4. **No Centralized Error Utilities**: Error display logic is in `utils::error::show_snafu_error()` but not consistently used

## Improvement Areas

### 1. Memory Overhead Reduction

**Problem**: Using `backtrace: Option<Backtrace>` in every error variant adds ~8 bytes per variant even when backtraces aren't enabled.

**Solutions**:

#### Option A: Enable Implicit Backtrace Feature (Recommended)
Enable snafu's `backtrace-feature` which uses `std::backtrace::Backtrace` implicitly:

```toml
# Cargo.toml
[dependencies]
snafu = { version = "0.8.9", features = ["backtrace"] }
```

Then remove explicit backtrace fields:

```rust
// Before
#[derive(Debug, Snafu)]
pub enum Error {
    StdIo {
        source: std::io::Error,
        backtrace: Option<Backtrace>,  // Remove this
    },
}

// After
#[derive(Debug, Snafu)]
pub enum Error {
    StdIo {
        source: std::io::Error,
        // Backtrace is implicit - only captured when std::backtrace is enabled
    },
}
```

**Benefits**:
- Reduces memory by ~8 bytes per error variant
- Backtraces only captured when `RUST_BACKTRACE=1` is set
- Cleaner error definitions
- More idiomatic snafu usage

#### Option B: Use `#[snafu(implicit)]` Selectively
For errors that rarely need backtraces, use `#[snafu(implicit)]`:

```rust
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(implicit)]  // No backtrace for this variant
    NotFound,
    
    StdIo {
        source: std::io::Error,
        // Has implicit backtrace
    },
}
```

### 2. Improve Error Context and Display Messages

**Problem**: Many errors lack informative context, making debugging difficult.

**Solutions**:

#### A. Add Custom Display Messages
All error variants should have meaningful display messages:

```rust
// Before
#[derive(Debug, Snafu)]
pub enum Error {
    IoError {
        source: std::io::Error,
    },
}

// After
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(display("Failed to perform I/O operation: {source}"))]
    IoError {
        source: std::io::Error,
    },
}
```

#### B. Add Context Fields Where Valuable
For file operations, include paths:

```rust
// Before
#[derive(Debug, Snafu)]
pub enum Error {
    StdIo {
        source: std::io::Error,
    },
}

// After - Following ibd.rs pattern
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(display("Failed to read file: {}", path.display()))]
    ReadFile {
        source: std::io::Error,
        path: std::path::PathBuf,
    },
    
    #[snafu(display("Failed to write file: {}", path.display()))]
    WriteFile {
        source: std::io::Error,
        path: std::path::PathBuf,
    },
}
```

For parsing errors, include the problematic value:

```rust
#[snafu(display("Failed to parse '{value}' as {type_name}"))]
ParseValue {
    source: ParseFloatError,
    value: String,
    type_name: String,
}
```

### 3. Improve Error Hierarchy and Organization

**Problem**: Flat error enums in large modules become unwieldy. Module boundaries unclear.

**Solutions**:

#### A. Use Error Modules for Large Components
Following the pattern in `utils.rs` with `utils::path::Error`:

```rust
// Before: Large flat enum
pub enum Error {
    IoError { ... },
    ParseError { ... },
    ValidationError { ... },
    // ... 20+ variants
}

// After: Organized into submodules
pub mod io {
    #[derive(Debug, Snafu)]
    pub enum Error {
        #[snafu(display("Failed to read {}", path.display()))]
        Read { source: std::io::Error, path: PathBuf },
        // ... other I/O errors
    }
}

pub mod parse {
    #[derive(Debug, Snafu)]
    pub enum Error {
        // ... parsing errors
    }
}

// Top-level error aggregates them
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Io { source: io::Error },
    
    #[snafu(transparent)]
    Parse { source: parse::Error },
}
```

#### B. Use Error Enums at Module Boundaries
Good examples already exist (`share.rs`, `genome.rs` with transparent errors)

### 4. Backtrace Display Improvements

**Current**: The `show_snafu_error()` function in `utils::error` does frame filtering.

**Improvements**:

#### A. Make Filtering Configurable
```rust
pub struct BacktraceConfig {
    pub filter_stdlib: bool,
    pub filter_crates: bool,
    pub max_frames: Option<usize>,
}

impl Default for BacktraceConfig {
    fn default() -> Self {
        Self {
            filter_stdlib: true,
            filter_crates: true,
            max_frames: Some(20),
        }
    }
}

pub fn show_snafu_error_with_config<E>(e: E, config: &BacktraceConfig)
where
    E: ErrorCompat + AsErrorSource,
{
    // Enhanced implementation
}
```

#### B. Add Environment Variable Control
```rust
// Check ISHARE_BACKTRACE_FILTER=full|minimal|none
let filter_mode = std::env::var("ISHARE_BACKTRACE_FILTER")
    .unwrap_or_else(|_| "minimal".to_string());
```

### 5. Consistent Error Construction Patterns

**Problem**: Inconsistent use of `.context()` vs direct error construction.

**Guidelines**:

#### A. Use `.context()` for External Errors
```rust
// Good
std::fs::read_to_string(path)
    .context(ReadFileSnafu { path: path.to_owned() })?;
```

#### B. Use `ensure!()` for Validation
```rust
// Good
ensure!(bp <= chrlen, MapBpOutOfRangeSnafu { bp, chrlen });
```

#### C. Use `.fail()` for Context-Only Errors
```rust
// Good
NotEnoughItemSnafu.fail()?;
```

## Implementation Phases

### Phase 1: Enable Implicit Backtraces (Low Risk)
**Estimated Impact**: Reduces memory overhead, improves consistency

1. Update `ishare-lib/Cargo.toml`:
   ```toml
   snafu = { version = "0.8.9", features = ["backtrace"] }
   ```

2. Remove explicit `backtrace: Option<Backtrace>` fields from all error enums

3. Remove manual backtrace imports: `use std::backtrace::Backtrace;`

4. Test that backtraces still work with `RUST_BACKTRACE=1`

**Files affected**: ~47 Rust files

### Phase 2: Enhance Context in High-Traffic Errors (Medium Risk)
**Estimated Impact**: Significantly improves debugging experience

1. Start with file I/O errors in:
   - `io.rs`
   - `gmap.rs`
   - `genome.rs`
   - `vcf.rs`

2. Add `path: PathBuf` context to all file operations

3. Add custom display messages with path information

**Files affected**: ~10-15 files

### Phase 3: Improve Error Display Messages (Low Risk)
**Estimated Impact**: Better error messages for end users

1. Add `#[snafu(display())]` attributes to all error variants without them

2. Ensure messages follow a consistent format:
   - Start with action: "Failed to...", "Unable to...", "Invalid..."
   - Include relevant context: values, paths, types
   - Keep messages concise but informative

**Files affected**: All 47 files with errors

### Phase 4: Organize Large Error Enums (Medium Risk)
**Estimated Impact**: Better code organization, maintainability

1. Identify modules with >10 error variants

2. Group related errors into submodules

3. Use transparent errors at module boundaries

**Files affected**: ~5-8 large modules

### Phase 5: Enhance Backtrace Display (Low Risk)
**Estimated Impact**: Better debugging experience

1. Add configuration to `show_snafu_error()`

2. Add environment variable control

3. Document usage in README

**Files affected**: `utils.rs`, documentation

## Memory Impact Analysis

### Current State
- Each `Option<Backtrace>` field: ~8 bytes when None, ~kb when Some
- 99 backtrace fields across codebase
- Average error size: ~24-40 bytes per variant

### After Phase 1
- Remove 99 explicit backtrace fields
- Memory savings: ~800 bytes per error instance (when no backtrace)
- Backtraces only allocated when `RUST_BACKTRACE=1` is set
- Average error size: ~16-32 bytes per variant

## Testing Strategy

### For Each Phase

1. **Compilation Test**: Ensure all code compiles after changes
2. **Unit Tests**: Run existing unit tests to catch regressions
3. **Integration Tests**: Run integration tests if available
4. **Manual Testing**: Test error scenarios with `RUST_BACKTRACE=1`
5. **Error Display Test**: Verify error messages are clear and informative

### Specific Tests

```rust
#[test]
fn test_error_has_backtrace_when_enabled() {
    std::env::set_var("RUST_BACKTRACE", "1");
    let err = some_error_causing_function();
    assert!(ErrorCompat::backtrace(&err).is_some());
}

#[test]
fn test_error_display_includes_context() {
    let err = ReadFileSnafu { 
        path: PathBuf::from("/test/path") 
    }.fail::<()>().unwrap_err();
    let msg = format!("{}", err);
    assert!(msg.contains("/test/path"));
}
```

## Best Practices Going Forward

1. **Always add display messages** for new error variants
2. **Include relevant context** (paths, values, types) in error structs
3. **Use transparent errors** for module boundaries
4. **Don't add explicit backtrace fields** - rely on implicit backtraces
5. **Test error messages** during code review
6. **Document error handling patterns** in CONTRIBUTING.md

## Examples

### Before and After: io.rs

```rust
// BEFORE
use std::backtrace::Backtrace;

#[derive(Snafu, Debug)]
#[snafu(visibility(pub(crate)))]
pub enum Error {
    Downcast {
        backtrace: Option<Backtrace>,
    },
    StdIo {
        source: std::io::Error,
        backtrace: Option<Backtrace>,
    },
    Parquet {
        source: parquet::errors::ParquetError,
        backtrace: Option<Backtrace>,
    },
}

// AFTER
#[derive(Snafu, Debug)]
#[snafu(visibility(pub(crate)))]
pub enum Error {
    #[snafu(display("Failed to downcast Arrow array to expected type"))]
    Downcast,
    
    #[snafu(display("I/O operation failed: {source}"))]
    StdIo {
        source: std::io::Error,
    },
    
    #[snafu(display("Parquet operation failed: {source}"))]
    Parquet {
        source: parquet::errors::ParquetError,
    },
}
```

### Before and After: gmap.rs with Rich Context

```rust
// BEFORE
#[derive(Debug, Snafu)]
pub enum Error {
    IoError {
        source: std::io::Error,
        backtrace: Option<Backtrace>,
    },
    NotEnoughItem {
        backtrace: Option<Backtrace>,
    },
}

// AFTER
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(display("Failed to read genetic map file: {}", path.display()))]
    ReadMapFile {
        source: std::io::Error,
        path: PathBuf,
    },
    
    #[snafu(display("Failed to write genetic map file: {}", path.display()))]
    WriteMapFile {
        source: std::io::Error,
        path: PathBuf,
    },
    
    #[snafu(display("Genetic map has insufficient data points"))]
    NotEnoughItem,
}
```

## Risks and Mitigation

### Risk 1: Breaking Changes
**Mitigation**: This is an internal change. External API remains the same. Test thoroughly.

### Risk 2: Performance Impact
**Mitigation**: Implicit backtraces have minimal overhead when disabled. Measure before/after.

### Risk 3: Existing Error Handling Code
**Mitigation**: Most code uses `.context()` which continues to work. Search for manual backtrace handling.

### Risk 4: Error Message Changes
**Mitigation**: If tests depend on exact error messages, update them. Consider semantic testing instead.

## Conclusion

This plan provides a structured approach to improve error handling:
- **Reduced memory overhead** through implicit backtraces
- **Better debugging** with rich context and clear messages
- **Improved organization** with modular error hierarchies
- **Enhanced user experience** with informative error messages

The phases can be implemented incrementally with testing at each step, minimizing risk while delivering continuous improvements.
