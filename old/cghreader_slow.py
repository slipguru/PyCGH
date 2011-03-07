#-*- coding: utf-8 -*-

import numpy as np

class CGHReader(object):

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
            types = acgh.readline().strip() # not used
            fields = acgh.readline().strip().split(delimiter)[1:]
            cols = np.arange(len(fields))+1

            self._data = np.genfromtxt(acgh,
                                       dtype=None,
                                       names=fields,
                                       usecols=cols,
                                       delimiter=delimiter,
                                       autostrip=True)
            # dictionary mapping
            self.FEATURES = dict((key, self._data[key]) for key in fields)

    def _extract_info(self, info_dict, acgh):
        types = acgh.readline().strip().split(self.delimiter)[1:]
        info = acgh.readline().strip().split(self.delimiter)[1:]

        data = acgh.readline().strip().split(self.delimiter)[1:]
        for i, t, d in zip(info, types, data):
            info_dict[i] = self.TYPE_MAP[t](d)

        self._skip_separator(acgh, check='*')

    def _skip_separator(self, acgh, check=None):
        separator = acgh.readline().strip()
        if not check is None:
            assert separator == check


# Tests -----------------------------------------------------------------------
from nose.tools import *

class TestCGHReader(object):
    def setup(self):
        self.reader = CGHReader('test.txt')

    def test_acces_FEPARAMS(self):
        value = self.reader.FEPARAMS['Protocol_Name']
        assert_equals('CGH-v4_95_Feb07 (Read Only)', value)

        value = self.reader.FEPARAMS['Scan_NumChannels']
        assert_equals(2, value)

        value = self.reader.FEPARAMS['Scan_MicronsPerPixelY']
        assert_equals(5.0, value)

        assert_raises(KeyError, self.reader.STATS.__getitem__, 'fake')

        assert_equals(34, len(self.reader.FEPARAMS))

    def test_access_STATS(self):
        value = self.reader.STATS['gDarkOffsetAverage']
        assert_equals(20.185, value)

        value = self.reader.STATS['rLocalBGInlierNum']
        assert_equals(45140, value)

        assert_raises(KeyError, self.reader.STATS.__getitem__, 'fake')

        assert_equals(143, len(self.reader.STATS))

    def test_access_FEATURES(self):
        value = self.reader.FEATURES['LogRatio']
        assert_equals(45214, len(value))

        assert_raises(KeyError, self.reader.FEATURES.__getitem__, 'fake')

        assert_equals(40, len(self.reader.FEATURES))
