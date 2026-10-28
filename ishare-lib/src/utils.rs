pub mod path {
    use crate::error::{IshareError, Result};
    use error_stack::*;

    use std::path::*;

    /// create path to a file from path prefix and an exension name
    ///
    /// This function will remove all extensions from the prefix filename and
    /// then add the the given extension name. It does not touch the file
    /// system; use `output_from_prefix` for output files.
    pub fn from_prefix(prefix: impl AsRef<Path>, suffix: &str) -> Result<PathBuf> {
        // remove all extensions from path
        let mut o = prefix.as_ref().to_path_buf();
        while o.extension().is_some() {
            o.set_extension("");
        }

        // add suffix to the path name
        o.set_extension(suffix);
        Ok(o)
    }

    /// same as `from_prefix`, but also create the parent folders of the
    /// returned path if they do not exist
    pub fn output_from_prefix(prefix: impl AsRef<Path>, suffix: &str) -> Result<PathBuf> {
        let o = from_prefix(prefix, suffix)?;
        if let Some(parent) = o.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .change_context(IshareError::Utils)
                    .attach_with(|| format!("cannot create folder {}", parent.display()))?;
            }
        }
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

    #[test]
    fn test_output_from_prefix_creates_parent() {
        let tmp = tempfile::tempdir().unwrap();
        let prefix = tmp.path().join("a/b/x.cdf");
        let p = output_from_prefix(&prefix, "abc").unwrap();
        assert_eq!(p, tmp.path().join("a/b/x.abc"));
        assert!(tmp.path().join("a/b").is_dir());

        // from_prefix itself does not create folders
        let p2 = from_prefix(tmp.path().join("c/d/y"), "abc").unwrap();
        assert_eq!(p2, tmp.path().join("c/d/y.abc"));
        assert!(!tmp.path().join("c").exists());

        // a bare file name has no parent folder to create
        assert_eq!(output_from_prefix("z.cdf", "abc").unwrap(), Path::new("z.abc"));
    }
}

pub mod error {
    // use regex::Regex;
    // use snafu::{AsErrorSource, Backtrace, ErrorCompat};

    // #[derive(Debug)]
    // struct Frame {
    //     func: String,
    //     file: String,
    //     line: u32,
    // }
    //     fn extract_frames(bt: &Backtrace) -> Vec<Frame> {
    //         let bt_str = format!("{bt:?}");
    //         let re = match Regex::new(r#"fn: "([^"]+)", file: "([^"]+)", line: (\d+)"#) {
    //             Ok(re) => re,
    //             Err(_) => return Vec::new(), // Return empty vec if regex fails
    //         };

    //         re.captures_iter(&bt_str)
    //             .filter_map(|cap| {
    //                 let func = cap.get(1)?.as_str().to_string();
    //                 let file = cap.get(2)?.as_str().to_string();
    //                 let line = cap.get(3)?.as_str().parse().ok()?;
    //                 Some(Frame { func, file, line })
    //             })
    //             .collect()
    //     }

    //     pub fn show_snafu_error<E>(e: E)
    //     where
    //         E: ErrorCompat + AsErrorSource,
    //     {
    //         for (ic, c) in ErrorCompat::iter_chain(&e).enumerate() {
    //             if ic == 0 {
    //                 eprintln!("ERROR");
    //             }
    //             eprintln!("{ic:>4}: {c}");
    //         }
    //         if let Some(bt) = ErrorCompat::backtrace(&e) {
    //             eprintln!("BACKTRACE");
    //             let mut iframe = 0;
    //             for frame in extract_frames(bt) {
    //                 let file = &frame.file;
    //                 let func = &frame.func;
    //                 if file.contains("/rustc/")
    //                     || file.contains("crates.io")
    //                     || file.contains("toolchains")
    //                     || func.contains("as snafu::IntoError")
    //                 {
    //                     continue;
    //                 }
    //                 eprintln!(
    //                     "{iframe:>4}: {}\n        {:}:{}",
    //                     frame.func, frame.file, frame.line
    //                 );
    //                 iframe += 1;
    //             }
    //         }
    //     }
}
