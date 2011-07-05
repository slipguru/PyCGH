#-*- coding: utf-8 -*-
import itertools as it

import numpy as np

from cghutils import ArrayCGH

# TO BE REMOVED!!!! ##########################################################
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

logger.addHandler(handler)
###############################################################################

# Utility functions -----------------------------------------------------------
            # Agilent -> (conversion, dtype)
TYPE_MAP = {'text': (unicode, unicode),
            'float': (float, float),
            'integer': (int, int),
            'boolean': (lambda x: bool(int(x)), bool)}

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
                print 'Skipping unreadable lines nro %d' % ln
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


# Conversions algoritm --------------------------------------------------------
def _split_mapping(location):
    try:
        chr, interval = location.split(':')
        start, end = (int(x) for x in interval.split('-'))
    except ValueError: # unmapped/control probe
        return (AgilentCGH.INVALID_INT,
                AgilentCGH.INVALID_INT,
                AgilentCGH.INVALID_INT, True)

    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr == 'X':
        return 23, start, end, False
    elif chr=='Y':
        return 24, start, end, False
    else:
        try:
            return int(chr), start, end, False
        except ValueError: # unplaceable probe eg chrUn
            return (AgilentCGH.INVALID_INT,
                    AgilentCGH.INVALID_INT,
                    AgilentCGH.INVALID_INT, True)

# Main Class -------------------------------------------------------------------
class AgilentCGH(object):#ArrayCGH):

    INVALID_INT = -9999
    INVALID_FLOAT = np.nan
    INVALID_STRING = 'N/A'
    QC_FLAGS = ('gIsFeatNonUnifOL', 'rIsFeatNonUnifOL',
                'gIsFeatPopnOL', 'rIsFeatPopnOL',
                'gIsBGNonUnifOL', 'rIsBGNonUnifOL',
                'gIsBGPopnOL', 'rIsBGPopnOL',
                'gIsSaturated', 'rIsSaturated')
    NQC_FLAGS = ('gIsWellAboveBG', 'rIsWellAboveBG')

    #def __init__(self, *args, **kwargs):
        #return super(AgilentCGH, self).__init__(*args, **kwargs)

    @staticmethod
    def load(path, delimiter='\t', test_channel='r',
             fill_missings=False, qc_masking=False, release=None):

        if not test_channel in ('r', 'g'):
            raise ValueError("test_channel must be 'r' (default) or 'g'")

        with open(path, 'r') as acgh:
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
        if not release is None:
            rel_data = np.genfromtxt(release, dtype=None, usecols=(4, 1, 2, 3),
                                     names=('id', 'chr', 'sb', 'eb'))
            rel_data['sb'] += 1 # ucsc starts from zero
            rel_map = dict((id, '%s:%d-%d' % (c, s, e))
                           for id, c, s, e in rel_data)

            # If an id is missing in the release,
            # the probe is marked as unmapped
            locations = (rel_map.get(id, 'unmapped')
                         for id in features['ProbeName'])
        else:
            locations = features['SystematicName']

        # Chromosome Position extraction (X=23, Y=24)
        # Note that the control probe will be automatically removed
        # during this step
        loc_buff = zip(*(_split_mapping(x) for x in locations))

        # Data extraction
        data = it.chain([features[k] for k in agilent_names], loc_buff[:-1])
        mask = np.array(loc_buff[-1])

        # Agilent: outliers
        if qc_masking:
            qc_mask = np.zeros_like(mask)
            for flag in AgilentCGH.QC_FLAGS:
                np.logical_or(qc_mask, features[flag], qc_mask)
            for flag in AgilentCGH.NQC_FLAGS:
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
                [AgilentCGH.INVALID_STRING]*missing_len,    #id
                list(missing_rows),
                list(missing_cols),
                [AgilentCGH.INVALID_FLOAT]*missing_len,     #reference_signal
                [AgilentCGH.INVALID_FLOAT]*missing_len,     #test_signal
                [AgilentCGH.INVALID_INT]*missing_len,       #chromosome
                [AgilentCGH.INVALID_INT]*missing_len,       #start_base
                [AgilentCGH.INVALID_INT]*missing_len,       #end_base
            )

            # Merging
            full_data = list()
            for d, md in zip(data, missing_data):
                full_data.append(np.r_[d, md])
            data = full_data
            mask = np.r_[mask, np.array([True]*missing_len)]

        # Creation and dinamyc attachment of useful informations
        aCGH = ArrayCGH(*data, mask=mask)
        aCGH.TEST_CHANNEL = test_channel
        aCGH.PARAMS = params
        aCGH.STATS = stats
        aCGH.FEATURES = features
        aCGH.NAMES_MAP = dict(zip(ArrayCGH.COL_NAMES[:6], agilent_names))

        return aCGH
