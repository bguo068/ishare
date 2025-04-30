# Developer Environment Setup

This document provides instructions for setting up a consistent and fully
functional development environment for the ishare project.

The project is written in Rust and relies on native libraries like `htslib` and its Rust binding `rust-htslib`, which
require a C compiler such as clang (which `bindgen` uses). Development can be done on:

- macOS
- Ubuntu Linux
- Windows (via WSL). [WSL setup guide](https://learn.microsoft.com/en-us/windows/wsl/install)

We recommend using VSCode with the Rust extension for ease in setup and the best development experience.

## 📦 Prerequisites

1. Install Rust (via rustup)

This installs the latest stable Rust compiler along with cargo, the Rust package manager:

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

2. Install Required Rust Components

These tools are essential for code completion, linting, and formatting:

```sh
rustup component add rust-analyzer
rustup component add clippy
rustup component add rustfmt
```

3. Install a C Compiler (Required for `rust-htslib`)

- macOS

```sh
xcode-select --install
```

- Ubuntu / WSL2

```sh
sudo apt update
sudo apt install -y build-essential pkg-config git clang
```

## 🧑‍💻 Setting Up VSCode (Recommended)

1. Install [Visual Studio Code](https://code.visualstudio.com/)
2. Install the Rust extension by the Rust Foundation (search [rust-analyzer](https://code.visualstudio.com/docs/languages/rust) in the Extensions tab)
3. Open the project folder and make sure you allow VSCode to install any recommended extensions

Recommended VSCode settings (.vscode/settings.json):

```json
{
  "rust-analyzer.check.command": "clippy",
  "editor.formatOnSave": true,
  "editor.codeActionsOnSave": {
    "source.fixAll": true
  }
}
```

## Fork and develop on your copy

Fork: go to https://github.com/bguo068/ishare.git and hit fork and choose an owner and a name to fork

Clone the forked repository and build:

```sh
# replace the url to that of your own fork
git clone https://github.com/bguo068/ishare.git
cd ishare
# pull submodules that contain some test data
git submodule update --init --recursive
```

## 🚀 Build and Run

```sh
# build for debuging
cargo build --bin gtencode
# build for release
cargo build --release --bin gtencode
```

## Write and run tests

Write tests: https://doc.rust-lang.org/book/ch11-00-testing.html

Run tests:

```sh
cargo test
```

## Compile and open documentation

```sh
cargo doc --open
```

This is popup a web page that contains documentation for all used crates as well as the current project.

## Push your changes

```
git add [FILES_YOUR_CHANGES]
git commit -m '[SHORT_MEANING_FULL_MESSAGE]'
git push
```

## Pull Request

once you push your changes to your fork.
Pull Request can made from your fork to the upstream repo: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork

### 📚 References

- Rust Language Book: https://doc.rust-lang.org/book/
- Cargo Book: https://doc.rust-lang.org/cargo/
- rust-htslib: https://github.com/rust-bio/rust-htslib
- arrow-rs: https://github.com/apache/arrow-rs
