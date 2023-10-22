pub mod path {
    use std::path::*;

    /// create path to a filefrom path prefix and an exension name
    ///
    /// This function will create all folders as needed and remove all
    /// extensions from the prefix filename and then add the the given extension
    /// name.
    pub fn from_prefix(
        prefix: impl AsRef<Path>,
        suffix: &str,
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        // if parental path does not exist, then try to create all folders needed
        if let Some(parent) = prefix.as_ref().parent() {
            if parent == Path::new("") {
                std::fs::create_dir_all(parent)?;
            }
        }
        // remove all extensions from path
        let mut o = prefix.as_ref().to_path_buf();
        while let Some(_) = o.extension() {
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
