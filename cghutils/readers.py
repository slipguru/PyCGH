#-*- coding: utf-8 -*-
import itertools as it

import numpy as np

from cghutils import ArrayCGH

# Utility functions -----------------------------------------------------------
            # Agilent -> (conversion, dtype)
TYPE_MAP = {'text': (unicode, unicode),
            'float': (float, float),
            'integer': (int, int),
            'boolean': (lambda x: bool(int(x)), bool)}

INVALID_INT = -9999
INVALID_FLOAT = np.nan
INVALID_STRING = 'N/A'

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

    for line in acgh:
        data = line.strip().split(delimiter)[1:]
        for i, t, d in zip(info, types, data):
            out.setdefault(i, []).append(TYPE_MAP[t][0](d))

    # Conversion in numpy array
    for i, t in zip(info, types):
        out[i] = np.asarray(out[i], dtype=TYPE_MAP[t][1])

    return out

def _return_headers(acgh, delimiter='\t'):
    types = acgh.readline().strip().split(delimiter)[1:]
    info = acgh.readline().strip().split(delimiter)[1:]
    return types, info


# Conversions algoritm --------------------------------------------------------
def _split_mapping(location):
    try:
        chr, interval = location.split(':')
        start, end = (int(x) for x in interval.split('-'))
    except ValueError: # unmapped/control probe
        return (INVALID_INT,INVALID_INT, INVALID_INT, True)

    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr == 'X':
        return 23, start, end, False
    elif chr=='Y':
        return 24, start, end, False
    else:
        return int(chr), start, end, False

# Main Classes ------------------------------------------------------------------
class AgilentCGH(ArrayCGH):

    def __init__(self, *args, **kwargs):
        return super(AgilentCGH, self).__init__(*args, **kwargs)

    @staticmethod
    def load(path, delimiter='\t', test_channel='r'):

        if not test_channel in ('r', 'g'):
            raise ValueError("test_channel must be 'r' (default) or 'g'")

        with open(path, 'r') as acgh:
            # Reading FEPARAMS
            params = _read_info_line(acgh, delimiter)
            # Reading STATS
            stats = _read_info_line(acgh, delimiter)
            # Reading FEATURES
            features = _read_info_block(acgh, delimiter)

        # Mapping between Agilent Names and aCGH.COL_NAMES
        agilent_names = ['ProbeName', 'Row', 'Col']
        if test_channel == 'r':
            agilent_names.extend(['gMedianSignal', 'rMedianSignal'])
        elif test_channel == 'g':
            agilent_names.extend(['rMedianSignal', 'gMedianSignal'])

        # Chromosome Position extraction (X=23, Y=24)
        locations = features['SystematicName']
        loc_buff = zip(*(_split_mapping(x) for x in locations))

        # Data extraction
        data = it.chain([features[k] for k in agilent_names], loc_buff[:-1])

        aCGH = AgilentCGH(data, mask=loc_buff[-1])

        # Dinamyc attachment of useful informations
        aCGH.TEST_CHANNEL = test_channel
        aCGH.PARAMS = params
        aCGH.STATS = stats
        aCGH.FEATURES = features
        aCGH.NAMES_MAP = dict(zip(ArrayCGH.COL_NAMES[:6], agilent_names))

        return aCGH


    def _extract_array(self):
        agilent_names = ['Row', 'Col', 'ProbeName',
                         'rMedianSignal', 'gMedianSignal',
                         'rBGMedianSignal', 'gBGMedianSignal',
                         'LogRatio']

        if self._test_channel == 'r':
            array_names =  ['row', 'col', 'id',
                            'test_signal', 'ref_signal',
                            'test_bg', 'ref_bg', 'ratio']
        elif self._test_channel == 'g':
            array_names =  ['row', 'col', 'id',
                            'ref_signal', 'test_signal',
                            'ref_bg', 'test_bg', 'ratio']

        # (chr, start, end) extracted from SystematicName
        # valid: not (unmapped , controls, OL, poor quality)
        array_names.extend(['chromosome', 'start_base', 'end_base', 'valid'])

        # Data lenght
        data_len = len(self._features['Row'])
        try:
            num_rows = self._params['Grid_NumRows']
            num_cols = self._params['Grid_NumCols']
        except KeyError:
            num_rows = self._features['Row'].max()
            num_cols = self._features['Col'].max()
        full_data_len = num_rows*num_cols

        # Chromosome Position extraction (X=23, Y=24)
        locations = self._features['SystematicName']
        loc_buff = it.izip(*(_split_mapping(x) for x in locations))

        # Data extraction
        data = it.chain([self._features[k] for k in agilent_names], loc_buff)
        rdata = np.rec.fromarrays(data, names=array_names).view(np.ndarray)

        # Agilent: outliers
        for flag in ('gIsFeatNonUnifOL', 'rIsFeatNonUnifOL',
                     'gIsFeatPopnOL', 'rIsFeatPopnOL',
                     'gIsBGNonUnifOL', 'rIsBGNonUnifOL',
                     'gIsBGPopnOL', 'rIsBGPopnOL'):
            # in-place update
            np.logical_and(rdata['valid'], ~self._features[flag], rdata['valid'])

        # Agilent: good quality probes
        np.logical_and(rdata['valid'],
                       np.logical_and(self._features['gIsWellAboveBG'],
                                      self._features['rIsWellAboveBG']),
                       rdata['valid'])

        # Missing data
        found_coords = set(tuple(x) for x in rdata[['row', 'col']])
        expexted_coords = set(it.product(xrange(1, num_rows+1), xrange(1, num_cols+1)))
        missing_rows, missing_cols = zip(*(expexted_coords - found_coords))

        missing_float = [INVALID_FLOAT]*len(missing_rows)
        missing_str = [INVALID_STRING]*len(missing_rows)
        missing_int = [INVALID_INT]*len(missing_rows)
        missing_bool = [False]*len(missing_rows)
        missing_data = [np.asarray(x) for x in (missing_rows,
                                                missing_cols,
                                                missing_str,    #id
                                                missing_float,  #ref_signal
                                                missing_float,  #sample_signal
                                                missing_float,  #ref_bg
                                                missing_float,  #sample_bg
                                                missing_float,  #ratio
                                                missing_int,    #chromosome
                                                missing_int,    #start_base
                                                missing_int,    #end_base
                                                missing_bool)]  #valid

        missing_rdata = np.rec.fromarrays(missing_data,
                                          names=array_names,
                                          dtype=rdata.dtype).view(np.ndarray)

        full_rdata = np.r_[rdata, missing_rdata]
        full_rdata.sort(order=('row', 'col')) # Default ordering

        return full_rdata
