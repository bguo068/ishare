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

### ⚠️ IMPORTANT: Phase 1 Has Been Revised

**See `PHASE1_DETAILED_ANALYSIS.md` for detailed testing and analysis.**

The original Phase 1 proposal to enable "implicit backtraces" by removing `backtrace: Option<Backtrace>` fields **will not work**. Testing confirmed that removing these fields eliminates backtrace capture entirely, even with `RUST_BACKTRACE=1` set.

### Phase 1A: Add Context and Display Messages (Low Risk) - REVISED
**Estimated Impact**: Significantly improves debugging without losing backtraces

1. Add `#[snafu(display())]` attributes to all error variants

2. Add context fields (paths, values, types) to errors where helpful:
   ```rust
   // Before
   IoError {
       source: std::io::Error,
       backtrace: Option<Backtrace>,
   }
   
   // After
   #[snafu(display("Failed to read file: {}", path.display()))]
   ReadFile {
       source: std::io::Error,
       path: PathBuf,
       backtrace: Option<Backtrace>,  // Keep this!
   }
   ```

3. **Keep all existing backtrace fields** - do not remove them

4. Test error messages are clear and informative

**Files affected**: All 47 files with error enums

### Phase 1B: Selective Backtrace Optimization (Optional, Medium Risk)
**Estimated Impact**: Small memory savings for specific error types

Only remove backtrace fields from:
- Simple validation errors (e.g., "value out of range")
- Errors that are self-explanatory with context alone
- High-frequency leaf errors where memory matters

**Requires case-by-case analysis and testing**

**Files affected**: ~10-15 files after careful analysis

### Phase 2: Organize Large Error Enums (Medium Risk)
**Estimated Impact**: Better code organization, maintainability

1. Identify modules with >10 error variants

2. Group related errors into submodules

3. Use transparent errors at module boundaries

**Files affected**: ~5-8 large modules

### Phase 3: Enhance Backtrace Display (Low Risk)
**Estimated Impact**: Better debugging experience

1. Add configuration to `show_snafu_error()`

2. Add environment variable control

3. Document usage in README

**Files affected**: `utils.rs`, documentation

## Memory Impact Analysis

### Current State
- Each `Option<Backtrace>` field: ~8-16 bytes when None, several KB when Some
- 99 backtrace fields across codebase
- Average error size: ~24-56 bytes per variant (depending on other fields)

### After Revised Phase 1A
- **No memory changes** - all backtrace fields retained
- Focus on better error context and messages
- Debugging capability fully preserved

### After Optional Phase 1B (Selective Removal)
- If ~30 simple errors have backtraces removed: ~240-480 bytes saved per error instance
- Complex/external errors keep backtraces
- Trade-off: small memory savings vs. selective loss of debug info

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

1. **Always add display messages** for new error variants using `#[snafu(display())]`
2. **Include relevant context** (paths, values, types) in error structs
3. **Use transparent errors** for module boundaries
4. **Keep `backtrace: Option<Backtrace>` fields** for most errors - the memory cost is worth the debugging capability
5. **Only remove backtraces** from simple, self-explanatory validation errors after careful consideration
6. **Test error messages** and backtrace capture during code review
7. **Document error handling patterns** in CONTRIBUTING.md

## Risks and Mitigation (Updated)

### Risk 1: Breaking Changes
**Mitigation**: Phase 1A only adds fields and attributes. External API remains the same. Test thoroughly.

### Risk 2: Loss of Debugging Information (CRITICAL)
**Previous Plan Risk**: Removing all backtraces would eliminate stack traces.
**Mitigation**: Revised plan keeps backtraces. Only consider selective removal in Phase 1B with testing.

### Risk 3: Existing Error Handling Code
**Mitigation**: Most code uses `.context()` which continues to work. No changes needed.

### Risk 4: Error Message Changes
**Mitigation**: Adding display messages improves error output. Update tests if needed.

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

// AFTER (Phase 1A - Add context and messages, KEEP backtraces)
use std::backtrace::Backtrace;

#[derive(Snafu, Debug)]
#[snafu(visibility(pub(crate)))]
pub enum Error {
    #[snafu(display("Failed to downcast Arrow array to expected type"))]
    Downcast {
        backtrace: Option<Backtrace>,  // Kept!
    },
    
    #[snafu(display("I/O operation failed: {source}"))]
    StdIo {
        source: std::io::Error,
        backtrace: Option<Backtrace>,  // Kept!
    },
    
    #[snafu(display("Parquet operation failed: {source}"))]
    Parquet {
        source: parquet::errors::ParquetError,
        backtrace: Option<Backtrace>,  // Kept!
    },
}
```

### Before and After: gmap.rs with Rich Context (Revised)

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

// AFTER (Phase 1A - Add context, keep backtraces)
use std::path::PathBuf;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(display("Failed to read genetic map file: {}", path.display()))]
    ReadMapFile {
        source: std::io::Error,
        path: PathBuf,
        backtrace: Option<Backtrace>,  // Kept!
    },
    
    #[snafu(display("Failed to write genetic map file: {}", path.display()))]
    WriteMapFile {
        source: std::io::Error,
        path: PathBuf,
        backtrace: Option<Backtrace>,  // Kept!
    },
    
    #[snafu(display("Genetic map has insufficient data points"))]
    NotEnoughItem {
        backtrace: Option<Backtrace>,  // Kept for now
    },
}

// OPTIONAL Phase 1B - Remove backtrace from simple validation error
#[derive(Debug, Snafu)]
pub enum Error {
    // ... file errors keep backtraces ...
    
    #[snafu(display("Genetic map has insufficient data points"))]
    NotEnoughItem {
        // Backtrace removed - validation error, context is clear
    },
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

## Conclusion (REVISED)

**The original plan has been significantly revised based on testing.**

### Key Changes:
- **Phase 1 corrected**: Cannot remove backtrace fields without losing backtraces
- **New focus**: Add context and display messages while keeping backtraces
- **Memory trade-off acknowledged**: Debugging capability > memory savings

### This Revised Plan Provides:
- **Better debugging** with rich context (paths, values) and clear messages ✅
- **Preserved stack traces** for all errors that currently have them ✅
- **Improved organization** with modular error hierarchies (future phases) ✅
- **Enhanced user experience** with informative error messages ✅
- **Optional memory optimization** for simple validation errors only (Phase 1B) ⚠️

### Implementation Strategy:
1. **Start with Phase 1A**: Add context and messages (safe, high value)
2. **Test thoroughly**: Ensure backtraces still work with `RUST_BACKTRACE=1`
3. **Consider Phase 1B**: Only after 1A is complete and tested
4. **Proceed incrementally**: Each phase with testing, minimizing risk

See `PHASE1_DETAILED_ANALYSIS.md` for detailed test results and evidence.
