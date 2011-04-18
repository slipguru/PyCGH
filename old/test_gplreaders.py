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
