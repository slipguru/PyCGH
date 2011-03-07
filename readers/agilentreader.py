#-*- coding: utf-8 -*-

import numpy as np

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
            self.FEPARAMS = dict()
            self._extract_info(self.FEPARAMS, acgh)
            # Reading STATS
            self.STATS = dict()
            self._extract_info(self.STATS, acgh)

            # Reading FEATURES
            self.FEATURES = dict()
            self._extract_info(self.FEATURES, acgh, while_eol=True)

    def _extract_info(self, info_dict, acgh, while_eol=False):
        types = acgh.readline().strip().split(self.delimiter)[1:]
        info = acgh.readline().strip().split(self.delimiter)[1:]

        if not while_eol:
            data = acgh.readline().strip().split(self.delimiter)[1:]
            for i, t, d in zip(info, types, data):
                info_dict[i] = self.TYPE_MAP[t](d)
            self._skip_separator(acgh, check='*')
        else:
            for line in acgh:
                data = line.strip().split(self.delimiter)[1:]
                for i, t, d in zip(info, types, data):
                    info_dict.setdefault(i, []).append(self.TYPE_MAP[t](d))

            # Conversion in numpy array
            for i, t in zip(info, types):
                conv = bool if t == 'boolean' else self.TYPE_MAP[t]
                info_dict[i] = np.asarray(info_dict[i], dtype=conv)

    def _skip_separator(self, acgh, check=None):
        separator = acgh.readline().strip()
        if not check is None:
            assert separator == check
