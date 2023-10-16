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
        ibd_formats: str,
    ) -> "PyIbdSet": ...
    def get_coverage(self, step: int) -> Dict[str, ndarray]:
        """
        get coverage

        return a python dictionary, which has keys: 'ChromosomeId', 'ChrPos',
        "GwPos", "Coverage" and each of the corresponding values is a 1-d numpy
        array. The returned dictionary can be convert to Pandas.DataFrame cheaply
        """
        ...
