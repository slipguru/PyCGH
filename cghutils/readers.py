#-*- coding: utf-8 -*-

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

        self.TYPE_MAP = {'text': unicode,
            'float': float,
            'integer': int,
            'boolean': lambda x: bool(int(x))}

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
