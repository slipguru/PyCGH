from nose.tools import *
import os

import numpy as np

from cghutils.readers import AgilentReader, GPLReader

class TestAgilentReader(object):
    def setup(self):
        par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = AgilentReader(os.path.join(par_dir, 'test_agilent.txt'))

    def test_acces_params(self):
        value = self.reader.param('Protocol_Name')
        assert_equals('CGH-v4_95_Feb07 (Read Only)', value)

        value = self.reader.param('Scan_NumChannels')
        assert_equals(2, value)

        value = self.reader.param('Scan_MicronsPerPixelY')
        assert_equals(5.0, value)

        assert_raises(KeyError, self.reader.param, 'fake')

    def test_params_keys(self):
        params = self.reader.params_list()
        assert_equals(34, len(params))

    def test_access_stats(self):
        value = self.reader.stat('gDarkOffsetAverage')
        assert_equals(20.185, value)

        value = self.reader.stat('rLocalBGInlierNum')
        assert_equals(45140, value)

        assert_raises(KeyError, self.reader.stat, 'fake')

    def test_stats_keys(self):
        stats = self.reader.stats_list()
        assert_equals(143, len(stats))

    def test_access_features(self):
        value = self.reader.feature('LogRatio')
        assert_equals(4460, len(value))

        assert_raises(KeyError, self.reader.feature, 'fake')

    def test_features_keys(self):
        features = self.reader.features_list()
        assert_equals(40, len(features))

    def test_toarray(self):
        value = self.reader.toarray()
        assert_equals(45220, len(value))
        assert_equals(10, len(value[0]))

        # Find present data
        present = np.logical_and(value['row'] == 19, value['col'] == 64)
        assert_equals('A_14_P102318', value['id'][present])

        # Find missing data
        missing = np.logical_and(value['row'] == 18, value['col'] == 64)
        assert_equals(18, value['row'][missing])
        assert_equals(64, value['col'][missing])
        assert_equals('N/A', value['id'][missing])
        assert_true(np.isnan(value['ref_signal'][missing]))
        assert_true(np.isnan(value['sample_bg'][missing]))
        assert_equals(-9999, value['chromosome'][missing])


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
