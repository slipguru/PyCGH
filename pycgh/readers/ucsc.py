from ..utils import _file_handle, location_normalization
from ..datatypes.arraycgh import ArrayCGH

# filter valid return only valid probes
def ucsc_mapping(ucsc_file, filter_valid=False):
    """Reads an input file describing a chip and returns a dictionary mapping the probes included in the specified chip to their relative base pairs.
    
    The input file must contain rows organized in 5 tab-separated columns, representing:
    
    #. Indexing field to speed chromosome range queries.
    #. Reference sequence chromosome or scaffold
    #. Start position in chromosome
    #. End position in chromosome
    #. Name of the probe
    
    Example: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/agilentCgh4x44k.txt.gz
    
    Parameters
    ----------
    
    ucsc_file : str
        The path to the file containing the chip specification.
        
    filter_valid : bool, optional (default: False)
        Whether to include or not lines with possibly missing values. If set to False (default), returns a list of 4-tuples (chromosome number, start base, end base, a flag indicating if the probe has been correctly read), if set to True in each tuple the last value is omitted.
    
    Returns
    -------
    
    ucsc_mapping : list
        The list of tuples representing the probes on the chip.
    
    """
    
    ucsc_mapping = dict()
    for line in _file_handle(ucsc_file):
        _, chr, sb, eb, id = line.split()
        value = location_normalization(chr, int(sb)+1, int(eb))
        if not filter_valid:
            ucsc_mapping[id] = value
        elif not value[3]: # Only If valid
            ucsc_mapping[id] = value[:3]

    return ucsc_mapping

