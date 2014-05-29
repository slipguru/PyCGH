import numpy as np

from .ucsc import ucsc_mapping as ucsc_reader
from ..datatypes.arraycgh import ArrayCGH
from ..utils import _file_handle, split_location

TYPE_MAP = {'IMAGE_ID': unicode,
            'GENE_EXPR_OPTION': unicode,
            'SEQ_ID': unicode,
            'PROBE_ID': unicode,
            'POSITION': int,
            'X': int,
            'Y': int,
            'MATCH_INDEX': int,
            'SEQ_URL': unicode,
            'PM': float,
            'MM': float}

def _read_all_data(fileobj, header, delimiter):
    out = dict()
    for line in fileobj:
        data = line.strip().split(delimiter)
        for h, d in zip(header, data):
            out.setdefault(h, []).append(TYPE_MAP[h](d))

    # Conversion in numpy array
    for k in out:
        out[k] = np.asarray(out[k], dtype=TYPE_MAP[k])

    return out

def _read_pair(file_path, delimiter):
    with _file_handle(file_path) as pair_file:
        meta_line = pair_file.readline()[1:].strip().split(delimiter)
        meta = dict((k, v) for (k, v) in (x.split('=') for x in meta_line))

        header = pair_file.readline().strip().split(delimiter)
        data = _read_all_data(pair_file, header, delimiter)

    return meta, data

def _probe_to_pos(id):
    chr, pos = id.split('FS')
    chr = lower()
    pos = int(pos)
    return split_location('%s:%s-%s' % (chr, str(pos), str(pos)))

# Main Function ---------------------------------------------------------------
def nimblegen(test_path, reference_path, delimiter='\t', ucsc_mapping=None):
    """
    Parse a text file produced by a Nimblegen platform.
    
    Parameters
    ----------
    
    test_path : file or str
        Either the path to the file containing the intensity values for the test signal or the file handler itself.
    
    reference_path : file or str
        Either the path to the file containing the intensity values for the reference signal or the file handler itself.
        
    delimiter : str, optional (default: ``'\t'``)
        The string used to separate columns in the data file.
    
    ucsc_mapping : str or dict, optional (default: ``None``)
        Either the dictionary returned by :func:`pycgh.readers.ucsc_mapping`, which describes the chip being simulated (which probes are present on it), or the path to the file containing the mapping.
    
    Returns
    -------
    
    aCGH : :class:`pycgh.datatypes.ArrayCGH`
        The object representing the Array CGH loaded from the input file.
    """

    test_meta, test_data = _read_pair(test_path, delimiter)
    ref_meta, ref_data = _read_pair(reference_path, delimiter)
    
    if not np.all(test_data['PROBE_ID'] == ref_data['PROBE_ID']):
        raise ValueError('probe id mismatch between pair files')

    # Mapping between probe and chromosomal position
    if not ucsc_mapping is None:
        INVALID_PROBE_ID_VALUE = (ArrayCGH.MISSING_INT, # chr
                                  ArrayCGH.MISSING_INT, # sb
                                  ArrayCGH.MISSING_INT, # eb
                                  True)                 # mask
    
        if str(ucsc_mapping) == ucsc_mapping:
            rel_map = ucsc_reader(ucsc_mapping)
        else:
            rel_map = ucsc_mapping
        locations = (rel_map.get(id, INVALID_PROBE_ID_VALUE)
                     for id in test_data['PROBE_ID'])
    else:
        # Standard mapping included into nimblegen file
        locations = (_probe_to_pos(x) for x in test_data['PROBE_ID'])
    
    chromosome, startchr, endchr, mask = zip(*locations)

    # Creation and dinamyc attachment of useful informations
    aCGH = ArrayCGH(id=test_data['PROBE_ID'],
                    row=test_data['Y'],
                    col=test_data['X'],
                    reference_signal=ref_data['PM'],
                    test_signal=test_data['PM'],
                    chromosome=chromosome,
                    start_base=startchr,
                    end_base=endchr,
                    mask=mask)
    aCGH.META_TEST = test_meta
    aCGH.META_REF = ref_meta
    #aCGH.NAMES_MAP = dict(zip(ArrayCGH.COL_NAMES[:6],
    #                          ('PROBE_ID', 'Y', 'X', 'PM', 'PM', 'PROBE_ID'))
    #                     )

    return aCGH



