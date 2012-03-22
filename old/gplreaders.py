#-*- coding: utf-8 -*-
import itertools as it

import numpy as np

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

###########à

from nose.tools import *
import os

import numpy as np

import cghutils.readers as cghr
from cghutils.readers import AgilentReader, GPLReader

class TestGPLReader(object):
    def setup(self):
        par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = GPLReader(os.path.join(par_dir, 'test_geogpl.txt'))

    def test_fields(self):
        fields = self.reader.fields_list()
        assert_equals(13, len(fields))

    def test_fields_description(self):
        value = self.reader.field_description('ID')
        assert_equals("Agilent feature number", value)

    def test_field_value(self):
        values = self.reader.field('ID')
        assert_equals(30, len(values))

        missing_values = self.reader.field('GENE_SYMBOL')
        assert_equals(30, len(missing_values))

    def test_fields_types(self):
        id_sum = self.reader.field('ID').sum()
        assert_equals(124159, id_sum)

        gene_symbols = self.reader.field('GENE_SYMBOL')
        assert_raises(TypeError, gene_symbols.sum)