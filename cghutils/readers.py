#-*- coding: utf-8 -*-
import itertools as it

import numpy as np

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


# Main Classes ------------------------------------------------------------------
class AgilentReader(object):

    def __init__(self, path, delimiter='\t'):

        self.delimiter = delimiter
        self.path = path

        with open(path, 'r') as acgh:
            # Reading FEPARAMS
            self._params = _read_info_line(acgh, self.delimiter)
            # Reading STATS
            self._stats = _read_info_line(acgh, self.delimiter)
            # Reading FEATURES
            self._features = _read_info_block(acgh, self.delimiter)

    def param(self, key):
        return self._params[key]

    def params_list(self):
        return self._params.keys()

    def stat(self, key):
        return self._stats[key]

    def stats_list(self):
        return self._stats.keys()

    def feature(self, key):
        return self._features[key]

    def features_list(self):
        return self._features.keys()

    def toarray(self):
        agilent_names = ['Row', 'Col', 'ProbeName',
                         'rMedianSignal', 'gMedianSignal',
                         'rBGMedianSignal', 'gBGMedianSignal']
        array_names = ['row', 'col', 'id',
                       'ref_signal', 'sample_signal',
                       'ref_bg', 'sample_bg',
                       'chromosome', 'start_base', 'end_base'] # extract from SystematicName

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
        loc_buff = [np.array([0]*data_len, dtype=int),
                    np.array([0]*data_len, dtype=int),
                    np.array([0]*data_len, dtype=int)]

        # Data extraction
        data = [self._features[k] for k in agilent_names] + loc_buff
        rdata = np.rec.fromarrays(data, names=array_names).view(np.ndarray)

        # Missing data
        import itertools as it
        found_coords = set(tuple(x) for x in rdata[['row', 'col']])
        expexted_coords = set(it.product(xrange(1, num_rows+1), xrange(1, num_cols+1)))
        missing_rows, missing_cols = zip(*(expexted_coords - found_coords))

        missing_float = [np.nan]*len(missing_rows)
        missing_str = ['N/A']*len(missing_rows)
        missing_int = [-9999]*len(missing_rows)
        missing_data = [np.asarray(x) for x in (missing_rows,
                                                missing_cols,
                                                missing_str,    #id
                                                missing_float,  #ref_signal
                                                missing_float,  #sample_signal
                                                missing_float,  #ref_bg
                                                missing_float,  #sample_bg
                                                missing_int,    #chromosome
                                                missing_int,    #start_base
                                                missing_int)]   #end_base

        missing_rdata = np.rec.fromarrays(missing_data,
                                          names=array_names,
                                          dtype=rdata.dtype).view(np.ndarray)

        full_rdata = np.r_[rdata, missing_rdata]
        full_rdata.sort(order=('row', 'col'))

        return full_rdata


class GPLReader(object):
    def __init__(self, path, delimiter='\t'):

        self.delimiter = delimiter
        self.path = path

        self._fields = dict()
        self._fields_values = dict()

        with open(path, 'r') as gplfile:
            # Reading headers descriptions
            for line in gplfile:
                if line.startswith('#'):
                    field, description = [x.strip() for x in line.split('=')]
                    self._fields[field[1:]] = description
                else:
                    header = [x.strip() for x in line.split(delimiter)]
                    assert list(sorted(header)) == list(sorted(self._fields.keys()))
                    break

            # Reading columns
            for line in gplfile:
                values = [x.strip() for x in line.split(delimiter)]

                if values == ['']: continue # skip empty lines

                # until shortest iterable is exausted
                for k, v in it.izip(header, values):
                    self._fields_values.setdefault(k, []).append(v)

        # Convertion in numpy objects
        # Try int, then float and then unicode
        for k in self._fields_values:
            try:
                out = np.asarray(self._fields_values[k], dtype=int)
            except ValueError:
                try:
                    out = np.asarray(self._fields_values[k], dtype=float)
                except ValueError:
                    out = np.asarray(self._fields_values[k], dtype=unicode)

            self._fields_values[k] = out

    def fields_list(self):
        return self._fields

    def field_description(self, key):
        return self._fields[key]

    def field(self, key):
        return self._fields_values[key]
