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
