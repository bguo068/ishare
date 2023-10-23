from typing import List, Tuple, Dict
from numpy import ndarray

class PyIbdSet:
    """
    A class representing IbdSet
    """

    def __init__(
        self,
        chromsize: List[int],
        chromnames: List[str],
        samples: List[str],
        gmaps: List[List[Tuple[int, float]]],
        ibd_files: List[str],
        ibd_format: str,
    ) -> "PyIbdSet": 
        """
        Create a PyIbdSet Object

        #### Arguments
        - `chromsize`: list of integers indicating the length of chromosomes
        - `chromnames`: list of strings indicating the name of chromosomes, in
        the same order as in chromsize.
        - `samples`: list of strings indicating the names of samples
        - `gmaps`: list of list of  of tuples indicating the genetic maps of
        chromosomes, in the same order as in chromsize. The element of the outer
        list is an inner list for a single chromosome. The element of inner list
        is a 2-tuple, (bp, cM), bp: base pair position, cM the corresonding cM
        coordinates corresponding to bp position
        - `ibd_files`: a list of IBD files. If `ibd_format` is tskibd, each IBD
        file must be a path to IBD file cooresponding to a chromosome (the same
        order as `chromsize`). As the chromosome id is not part of the IBD file,
        and the order of the files indicate chromosome id. For other formats,
        such as `hapibd`, chromosme information is encoded in the IBD files, the
        IBD_files can be a list of any number of IBD files. They are just
        concatenated internally. 
        - `ibd_format`: can be one of the following `tskibd`, `hapibd` and `hmmibd`

        #### Return
        an `PyIbdSet` object


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

    def get_xirs(self, vcf_files:List[str], min_maf:float, method: int)-> Dict[str, ndarray]:
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
        cheaply
        """
        ...