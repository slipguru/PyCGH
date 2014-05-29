import itertools as it

import numpy as np

from .ucsc import ucsc_mapping as ucsc_reader
from ..datatypes.arraycgh import ArrayCGH
from ..utils import _file_handle, split_location


# Constants -------------------------------------------------------------------
QC_FLAGS = ('gIsFeatNonUnifOL', 'rIsFeatNonUnifOL',
            'gIsFeatPopnOL', 'rIsFeatPopnOL',
            'gIsBGNonUnifOL', 'rIsBGNonUnifOL',
            'gIsBGPopnOL', 'rIsBGPopnOL',
            'gIsSaturated', 'rIsSaturated')
"""
The list of fields...
"""

NQC_FLAGS = ('gIsWellAboveBG', 'rIsWellAboveBG')

            # Agilent -> (conversion, dtype)
TYPE_MAP = {'text': (unicode, unicode),
            'float': (float, float),
            'integer': (int, int),
            'boolean': (lambda x: bool(int(x)), bool)}

# Main Function ---------------------------------------------------------------
def agilent(path, delimiter='\t', test_channel='r',
            fill_missings=False, qc_masking=False, ucsc_mapping=None):
    """
    Parse a text file produced by an Agilent Technologies platform.
    
    Parameters
    ----------
    
    path : str
        The path to the file which stores Array CGH data.
    
    delimiter : str, optional (default: ``'\t'``)
        The string used to separate columns in the data file.
    
    test_channel : str, optional (default: ``'r'``)
        Which channel is associated to the test values, can be either ``'r'`` (default value, stands for *red*) or ``'g'`` (for *green*), depending on which fluorescent dye is used to highlight the two channels.
    
    fill_missings : bool, optional (default: ``False``)
    
    qc_masking : bool, optional (default: ``False``)
        If set to true, mask all values which for some reason do not pass a quality check, meaning that for that value at least one of the flags listed in constant ``QC_FLAGS`` is ``True`` OR at least one of the flags listed in constant ``NQC_FLAGS`` is ``False``.
    
    ucsc_mapping : str or dict, optional (default: ``None``)
        Either the dictionary returned by :func:`pycgh.readers.ucsc_mapping`, which describes the chip being simulated (which probes are present on it), or the path to the file containing the mapping.
    
    Returns
    -------
    
    aCGH : :class:`pycgh.datatypes.ArrayCGH`
        The object representing the Array CGH loaded from the input file.
    """

    if not test_channel in ('r', 'g'):
        raise ValueError("test_channel must be 'r' (default) or 'g'")

    with _file_handle(path) as acgh:
        # Reading FEPARAMS
        params = _read_info_line(acgh, delimiter)
        # Reading STATS
        stats = _read_info_line(acgh, delimiter)
        # Reading FEATURES
        features = _read_info_block(acgh, delimiter)

    # Data lenght
    data_len = len(features['Row'])
    try:
        num_rows = params['Grid_NumRows']
        num_cols = params['Grid_NumCols']
    except KeyError:
        num_rows = features['Row'].max()
        num_cols = features['Col'].max()

    # Mapping between Agilent Names and aCGH.COL_NAMES
    agilent_names = ['ProbeName', 'Row', 'Col']
    if test_channel == 'r':
        agilent_names.extend(['gMedianSignal', 'rMedianSignal'])
    elif test_channel == 'g':
        agilent_names.extend(['rMedianSignal', 'gMedianSignal'])

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
                     for id in features['ProbeName'])
    else:
        # Standard mapping included into agilent file
        locations = (split_location(x) for x in features['SystematicName'])

    # Chromosome Position extraction (X=23, Y=24)
    # Note that the control probe will be automatically removed
    # during this step
    loc_buff = zip(*locations)

    # Data extraction
    data = it.chain([features[k] for k in agilent_names], loc_buff[:-1])
    mask = np.array(loc_buff[-1])

    # Agilent: outliers
    if qc_masking:
        qc_mask = np.zeros_like(mask)
        for flag in QC_FLAGS:
            np.logical_or(qc_mask, features[flag], qc_mask)
        for flag in NQC_FLAGS:
            np.logical_or(qc_mask, ~features[flag], qc_mask)
        np.logical_or(mask, qc_mask, mask)

    # Fill missing data
    if fill_missings:
        # Missing data
        found_coords = set(zip(features['Row'], features['Col']))
        expected_coords = set(it.product(xrange(1, num_rows+1), xrange(1, num_cols+1)))
        missing_rows, missing_cols = zip(*(expected_coords - found_coords))
        missing_len = len(missing_rows)

        missing_data = (
            [ArrayCGH.MISSING_STRING]*missing_len,    #id
            list(missing_rows),
            list(missing_cols),
            [ArrayCGH.MISSING_FLOAT]*missing_len,     #reference_signal
            [ArrayCGH.MISSING_FLOAT]*missing_len,     #test_signal
            [ArrayCGH.MISSING_INT]*missing_len,       #chromosome
            [ArrayCGH.MISSING_INT]*missing_len,       #start_base
            [ArrayCGH.MISSING_INT]*missing_len,       #end_base
        )

        # Merging
        full_data = list()
        for d, md in zip(data, missing_data):
            full_data.append(np.r_[d, md])
        data = full_data
        mask = np.r_[mask, np.array([True]*missing_len)]

    # Creation and dinamyc linking of useful informations
    aCGH = ArrayCGH(*data, mask=mask)
    aCGH.TEST_CHANNEL = test_channel
    aCGH.PARAMS = params
    aCGH.STATS = stats
    aCGH.FEATURES = features
    aCGH.NAMES_MAP = dict(zip(ArrayCGH.COL_NAMES[:6], agilent_names))

    return aCGH

# Utility functions -----------------------------------------------------------
def _read_info_line(acgh, delimiter='\t'):
    out = dict()
    types, info = _return_headers(acgh, delimiter)

    data = acgh.readline().strip().split(delimiter)[1:]
    for i, t, d in zip(info, types, data):
        out[i] = TYPE_MAP[t][0](d)
    separator = acgh.readline().strip()
    assert separator == '*'

    return out

def _read_info_block(acgh, delimiter='\t'):
    out = dict()
    types, info = _return_headers(acgh, delimiter)

    num = 0
    for ln, line in enumerate(acgh):
        data = line.strip().split(delimiter)[1:]
        for i, t, d in zip(info, types, data):
            try:
                out.setdefault(i, []).append(TYPE_MAP[t][0](d))
            except Exception, e:
                print 'Skipping unreadable lines #%d' % ln
                for k in out:
                    if len(out[k]) == (num+1): #rollback
                        out[k] = out[k][:-1]
                break
        else:
            num += 1 #number of lines

    # Conversion in numpy array
    for i, t in zip(info, types):
        out[i] = np.asarray(out[i], dtype=TYPE_MAP[t][1])

    return out

def _return_headers(acgh, delimiter='\t'):
    types = acgh.readline().strip().split(delimiter)[1:]
    info = acgh.readline().strip().split(delimiter)[1:]
    return types, info
