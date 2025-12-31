use snafu::prelude::*;
use std::backtrace::Backtrace;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Path {
        // leaf
        source: path::Error,
        backtrace: Box<Option<Backtrace>>,
    },
}

pub mod path {

    use snafu::prelude::*;
    type Result<T> = std::result::Result<T, Error>;

    #[derive(Debug, Snafu)]
    pub enum Error {
        StdIo {
            source: std::io::Error,
            backtrace: Box<Option<Backtrace>>,
        },
    }

    use std::{backtrace::Backtrace, path::*};

    use snafu::Snafu;

    /// create path to a filefrom path prefix and an exension name
    ///
    /// This function will create all folders as needed and remove all
    /// extensions from the prefix filename and then add the the given extension
    /// name.
    pub fn from_prefix(prefix: impl AsRef<Path>, suffix: &str) -> Result<PathBuf> {
        // if parental path does not exist, then try to create all folders needed
        if let Some(parent) = prefix.as_ref().parent() {
            if parent == Path::new("") {
                std::fs::create_dir_all(parent).context(StdIoSnafu {})?;
            }
        }
        // remove all extensions from path
        let mut o = prefix.as_ref().to_path_buf();
        while o.extension().is_some() {
            o.set_extension("");
        }

        // add suffix to the path name
        o.set_extension(suffix);
        Ok(o)
    }

    #[test]
    fn test_from_prefix() {
        assert_eq!(
            Path::new("/tmp/1/2/x.abc"),
            from_prefix("/tmp/1/2/x.cdf", "abc").unwrap()
        );
        assert_eq!(
            Path::new("/tmp/1/2/x.abc"),
            from_prefix("/tmp/1/2/x.cdf.efg", "abc").unwrap()
        );
        assert_eq!(
            Path::new("/tmp/1/2/x.abc"),
            from_prefix("/tmp/1/2/x", "abc").unwrap()
        );
    }
}

pub mod error {
    use regex::Regex;
    use snafu::{AsErrorSource, Backtrace, ErrorCompat};

    #[derive(Debug)]
    struct Frame {
        func: String,
        file: String,
        line: u32,
    }
    fn extract_frames(bt: &Backtrace) -> Vec<Frame> {
        let bt_str = format!("{bt:?}");
        let re = match Regex::new(r#"fn: "([^"]+)", file: "([^"]+)", line: (\d+)"#) {
            Ok(re) => re,
            Err(_) => return Vec::new(), // Return empty vec if regex fails
        };

        re.captures_iter(&bt_str)
            .filter_map(|cap| {
                let func = cap.get(1)?.as_str().to_string();
                let file = cap.get(2)?.as_str().to_string();
                let line = cap.get(3)?.as_str().parse().ok()?;
                Some(Frame { func, file, line })
            })
            .collect()
    }

    pub fn show_snafu_error<E>(e: E)
    where
        E: ErrorCompat + AsErrorSource,
    {
        for (ic, c) in ErrorCompat::iter_chain(&e).enumerate() {
            if ic == 0 {
                eprintln!("ERROR");
            }
            eprintln!("{ic:>4}: {c}");
        }
        if let Some(bt) = ErrorCompat::backtrace(&e) {
            eprintln!("BACKTRACE");
            let mut iframe = 0;
            for frame in extract_frames(bt) {
                let file = &frame.file;
                let func = &frame.func;
                if file.contains("/rustc/")
                    || file.contains("crates.io")
                    || file.contains("toolchains")
                    || func.contains("as snafu::IntoError")
                {
                    continue;
                }
                eprintln!(
                    "{iframe:>4}: {}\n        {:}:{}",
                    frame.func, frame.file, frame.line
                );
                iframe += 1;
            }
        }
    }
}
