# `ishare`

`ishare` is a Rust crate designed to facilitate the analysis of rare-variant
sharing and identity-by-descent (IBD) sharing.

Currently, it's in pre-alpha stage.

## Rare Variant Sharing Analysis

This crate introduces a data structure and its corresponding algorithm for the
tabular encoding of rare variant genotype data. The tabular encoding leverages
the sparsity of minor alleles at sites with low minor allele frequencies,
effectively transforming a large site-oriented genotype matrix into a more
compact, manipulable, and accessible data structure. As a result, this design
enhances disk IO speeds, facilitates in-memory access to rare variants, and
permits large-scale rare variant sharing analysis.

- **Command Line Tool:** `gtencode`
- **Python Package (work in progress):** `isharepy`

## Identity-by-Descent Sharing Analysis

Included in this crate are the data structure and algorithms associated with IBD
sharing analysis. It serves as a reimplementation of functionalities present in
the Python package [`ibdutils`](https://github.com/bguo068/ibdutils) and the C++
library [`ibdtools`](https://github.com/umb-oconnorgroup/ibdtools). Key updates
from the previous implementations include:

- Use of genome-wide base-pair coordinates for IBD segment start and end fields,
as opposed to chromosomal positions. This facilitates easier IBD segment
manipulation across chromosomes.
- The [`itervaltree`](https://github.com/main--/rust-intervaltree) crate has
been integrated to improve various IBD processing and analysis functionalities,
enhancing code maintainability and readability.

- **Command Line Tool:** `ibdutils`
- **Python Package (work in progress):** `isharepy` for IBD-related classes.

## Ancestry-Specific IBD 

We've introduced a command line tool that integrates local ancestry information
from [`rfmix` v2](https://github.com/slowkoni/rfmix) with IBD segments to derive
ancestry-specific IBD. Refer to the binary `asibd` for this feature.

## Requirements

1. **Operating Systems:** Supported on Linux and MacOS. Windows is currently not
supported due to compatibility issues with the `rust-htslib` dependency.
2. **C Compiler:** This crate relies on bindings to C libraries.
3. **Rust:** Ensure `cargo` and the necessary toolchain are installed to compile
this crate and its associated binary executables.

## Compilation

```sh
cargo build --release --bin gtencode
cargo build --release --bin ibdutils
cargo build --release --bin asibd
```

Post-compilation, binaries are located in the target/release/ directory.

## Install the python package
1. make sure `cargo`, the necessary toolchain, c compilerare, python 3 are installed
2. use `pip` to install `maturin`, `numpy` and `pyarrow`
```
pip install maturin numpy pyarrow
```
3. install `isharepy` 
```
cd ishare_py
maturin develop --release
cd ..
```

# Usage
Instructions for each command line tool are as follows:

- `toolname --help`: Lists available subcommands. E.g., `gtencode --help`.
- `toolname subcommand --help`: Displays usage for a specific subcommand. 
E.g., `gtencode encode --help`.

# Test

Running tests requires cloning test data from `https://github.com/bguo068/testdata`, 
which has been added as  a submodule of `ishare` and thus can downloaded by running: 
`git submodule update --init --recursive`
