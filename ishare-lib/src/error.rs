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
