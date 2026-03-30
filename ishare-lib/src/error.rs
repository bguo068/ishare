use error_stack::Report;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum IshareError {
    #[error("dense genotype error")]
    CommonGenotype,
    #[error("rare genotype error")]
    RareGenotype,
    #[error("ibd error")]
    Ibd,
    #[error("ancestry-specific ibd error")]
    AsIbd,
    #[error("rv share error")]
    RvShare,
    #[error("genome error")]
    Genome,
    #[error("gmap error")]
    Gmap,
    #[error("individual error")]
    Individual,
    #[error("io error")]
    Io,
    #[error("site error")]
    Site,
    #[error("vcf error")]
    Vcf,
    #[error("unexpected empty Option value")]
    EmptyOption,
    #[error("runtime check error")]
    RuntimeCheck,
    #[error("utils error")]
    Utils,
    #[error("stats error")]
    Stats,
    #[error("matrix error")]
    Matrix,
    #[error("container error")]
    Container,
}

pub type Result<T> = std::result::Result<T, Report<IshareError>>;

#[derive(Debug)]
pub struct PathAttachment(std::path::PathBuf);

impl std::fmt::Display for PathAttachment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error related path: {}", self.0.display())
    }
}

impl From<&std::path::Path> for PathAttachment {
    fn from(path: &std::path::Path) -> Self {
        PathAttachment(path.into())
    }
}

impl From<std::path::PathBuf> for PathAttachment {
    fn from(path: std::path::PathBuf) -> Self {
        PathAttachment(path)
    }
}

impl From<String> for PathAttachment {
    fn from(path: String) -> Self {
        PathAttachment(std::path::PathBuf::from(path))
    }
}

impl From<&str> for PathAttachment {
    fn from(path: &str) -> Self {
        PathAttachment(std::path::PathBuf::from(path))
    }
}

impl From<&std::path::PathBuf> for PathAttachment {
    fn from(path: &std::path::PathBuf) -> Self {
        PathAttachment(path.clone())
    }
}
