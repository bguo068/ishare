from typing import List, Tuple, Dict
from numpy import ndarray

class GInfo:
    """
    A class representing GenomeInfo.
    """

    def __init__(
        name: str,
        chromsizes: List[int],
        chromnames: List[str],
        addition_chromname_2_id_map: Dict[str, int],
    ) -> "GInfo":
        """
        Create a GInfo Object

        #### Arguments
        - `chromsize`: list of integers indicating the length of chromosomes
        - `chromnames`: list of strings indicating the name of chromosomes, in
          the same order as in chromsize.
        - `addition_chromname_2_id_map`: a dictionary with key as str and value
          as int. This is can be useful. if you have differnt convention of
          chrnames that are mapped to this chrom idx. For example, you have
          chromnames = ["1", "2"], but also want to map "chr1" to index 0 and
          map "chr2" to index 1. This can be done by specifying
          `addition_chromname_2_id_map = {"chr1": 0, "chr2": 1}`
        #### Return
        an `GInfo` object
        """
        ...

class GMap:
    """
    A class representing GenomeInfo.
    """

    @staticmethod
    def from_list(lst: List[Tuple[int, int]], chrlen: int) -> "GMap":
        """
        Create a GMap for a single chromosome from a list of (bp, cM) pairs.

        #### Arguments
        - `lst`: a list a 2-tuple (bp, cM). Note bp should be 0-based.
        - `chrlen`: length of chromosome in bp.
        #### Return
        an `GMap` object
        """
        ...
    @staticmethod
    def from_plink_map_file(path: str, chrlen: int) -> "GMap":
        """
        Create a GMap for a single chromosome from reading a single-chromosome plink file

        #### Arguments
        - `path`: str, path to a single chromosome plink map file
        - `chrlen`: length of chromosome in bp.
        #### Return
        an `GMap` object
        """
        ...
    @staticmethod
    def from_constant_rate(rate: float, chrlen: int) -> "GMap":
        """
        Create a GMap for a single chromosome using constant rate

        #### Arguments
        - `rate`: float, rate of recombination in unit of "cM per bp"
        - `chrlen`: length of chromosome in bp.
        #### Return
        an `GMap` object
        """
        ...
    @staticmethod
    def from_list_of_gmap(lst: List[GMap]) -> "GMap":
        """
        Create a gneome-wide GMap by combining a list of a single chromosome GMap objects

        #### Arguments
        - `lst`: a list GMap object
        #### Return
        an `GMap` object
        """
        ...

class IBD:
    """
    A class representing IbdSet
    """

    def __init__(
        self,
        ginfo: GInfo,
        gmaps: GMap,
        ibd_files: List[str],
        ibd_format: str,
        samples: List[str],
    ) -> "IBD":
        """
        Create a IBD Object

        #### Arguments
        - `ginfo`: GInfo object
        - `gmaps`: GMap object
        - `ibd_files`: a list of IBD files. If `ibd_format` is tskibd, each IBD
        file must be a path to IBD file cooresponding to a chromosome (the same
        order as `chromsize`). As the chromosome id is not part of the IBD file,
        and the order of the files indicate chromosome id. For other formats,
        such as `hapibd`, chromosme information is encoded in the IBD files, the
        IBD_files can be a list of any number of IBD files. They are just
        concatenated internally.
        - `ibd_format`: can be one of the following `tskibd`, `hapibd` and `hmmibd`
        - `samples`: list of strings indicating the names of samples

        #### Return
        an `IBD` object
        """
        ...
    def get_coverage(self, step: int) -> Dict[str, ndarray]:
        """
        get coverage of the given IBD set

        #### Arguments
        - `step`: step size of the sampling point for coverage calculation

        #### Return
        a python dictionary, which has keys: 'ChromosomeId', 'ChrPos',
        "GwPos", "Coverage" and each of the corresponding values is a 1-d numpy
        array. The returned dictionary can be convert to Pandas.DataFrame cheaply
        """
        ...
    def get_xirs(
        self, vcf_files: List[str], min_maf: float, method: int
    ) -> Dict[str, ndarray]:
        """
        get Xir,S statistics table.

        #### Arguments
        - `vcf_files`: a list of vcf files, the length can be 1 or multiple. The
        purpose is to calculate of allele frequency from genotype data. The
        position of the sites will be canctenated from all VCF files. Only
        biallelic SNPs are considered.
        - `min_maf`: sites with minor allele freqeuncy lower than `min_maf` will
        be excluded
        - `method`: can be 1 or 2. `1` indicate the methods used in isorelate
        paper; `2` is an modified, experimental statistics.


        #### Return
        a python dictionary, which has keys: 'ChrId', 'ChrPos', "GwPos", "Raw",
        "Xirs", and "Pvalue", and each of the corresponding values is a 1-d
        numpy array. The returned dictionary can be convert to Pandas.DataFrame
        without much memory allocation.
        """
        ...

class RVar:
    """
    A class representing Rare variant genotype
    """

    @staticmethod
    def from_vcf(
        ginfo: GInfo,
        samples: List[str],
        vcf: str,
        chunksize: int,
        max_maf: float,
    ) -> "RVar":
        """
        Encode VCF/BCF file into tabular encoding format (REC) in RAM.

        #### Arguments
        - `ginfo`: GInfo object.
        - `samples`: a list of samples to work on. if the list is empty
        (length = 0), all samples in the vcf/BCF format will be used.
        - `vcf`: path to a VCF and BCF file
        - `chunksize`: window size in base pair for processing large vcf/bcf
        file in parallel.
        - `min_maf`: sites with minor allele freqeuncy lower than `min_maf` will
        be excluded

        #### Return
        return an `RVar` object
        """
        ...
    def __init__(self) -> "RVar":
        """
        Create an empty `RVar` file.

        Once created you can you the follow methods for reading different parts
        of encoded information:
        - `read_genotype`: for read genotype records from the `.rec` file
        - `read_sites`: read site information (minor and major alleles, site
        location) from the `.sit` file
        - `read_individuals`: read sample name information from `.ind` file
        """
        ...
    def read_genotype(self) -> None:
        """
        read genotype records into an existing `RVar` object. Any record currently
        held in the object will be cleared up before loading genotype records
        from the encoded file
        """
    def read_sites(self) -> None:
        """
        read sites info into an existing `RVar` object."""
    def read_inds(self) -> None:
        """
        read individual info into an existing `RVar` object."""
    @[staticmethod]
    def from_files(
        prefix: str, read_genotype: bool, read_sites: bool, read_individuals: bool
    ) -> "RVar":
        """
        read encoded files into a `RVar` object (without creating an empty RVar
        first)

        #### Arguments
        - `prefix`: prefix to the encoded file (any character after `.` in
        filename parts will be trimmed)
        - `read_genotype`: bool, indicating whether genotype records will be loaded.
        - `read_sites`: bool, indicating whether sites info will be loaded.
        - `read_individuals`: bool, indicating whether individuals' info will be loaded.
        #### Return
        return an `RVar` object
        """
        ...
    def get_genotype_table(self) -> Dict:
        """
        Convert `RVar` genotype into a 3-column table as a dictionary
        with keys `GenomeID`, `GenomeWidePosition` and `AlleleID`,
        (which can be easily converted to pd.DataFrame object)
        """
        ...
    def get_sites_table(self) -> Dict:
        """
        Convert `RVar` sites info into a 2-column table as a dictionary
        with keys `GenomeWidePosition` and `Alleles`,
        (which can be easily converted to pd.DataFrame object)
        """
        ...
    def calculate_jaccard_matrix(
        self,
        ids1: List[int] = None,
        ids2: List[int] = None,
        chunksize: int = None,
        copy_chunk: bool = None,
    ) -> Dict:
        """
        calculate genome-pair jaccard index

        #### Arguments
        - ids1: List, optional list of sample idx list.
        - ids2: List, optional list of sample idx list.
        - `chunksize`: int, indicate how many genomes to group together for
        parallelization.
        - `copy_chunk`: bool, indicating whether to make a copy of selected
        genome before jaccard matrix so genotype records for the selected genomes
        are sequential in RAM.
        #### Return
        return a dictionary with threes keys:
        - `RowIds`: 1-d ndarray indicating the genome id for each row
        - `ColumnIds`: 1-d ndarray indicating the genome id for each column
        - `Matrix`: 2-d ndarray array, the actually jaccard index result matrix
        """
        ...
