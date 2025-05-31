use std::{backtrace::Backtrace, num::ParseIntError};

use snafu::{ResultExt, Snafu};

#[derive(Debug, Snafu)]
enum MyError {
    IoError {
        source: std::io::Error,
        backtrace: Option<Backtrace>,
    },
    ParseInError {
        source: ParseIntError,
        backtrace: Option<Backtrace>,
    },
}

fn read_string_from_file(p: &str) -> Result<String, MyError> {
    std::fs::read_to_string(p).context(IoSnafu {})
}

fn parse_int(s: &str) -> Result<u32, MyError> {
    s.parse().context(ParseInSnafu {})
}

fn main() -> Result<(), MyError> {
    use ishare::utils::error::*;

    if let Err(e) = read_string_from_file("hello.txt") {
        show_snafu_error(e);
    }

    if let Err(e) = parse_int("123e") {
        show_snafu_error(e);
    }

    Ok(())
}
