from nose.tools import *
import os

import numpy as np

import cghutils.readers as cghr
from cghutils.readers import AgilentReader, GPLReader

class TestAgilentReader(object):
    def setup(self):
        self.par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = AgilentReader(os.path.join(self.par_dir, 'test_agilent.txt'))

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
        assert_equals(11, len(value[0]))

        # Find present data
        present = np.logical_and(value['row'] == 19, value['col'] == 64)
        assert_equals(19, value['row'][present])
        assert_equals(64, value['col'][present])
        assert_equals('A_14_P102318', value['id'][present])
        assert_equals(19, value['chromosome'][present])
        assert_equals(11455241, value['start_base'][present])
        assert_equals(11455291, value['end_base'][present])
        assert_true(~np.isnan(value['ref_signal'][present]))
        assert_true(~np.isnan(value['test_bg'][present]))

        # Find missing data
        missing = np.logical_and(value['row'] == 18, value['col'] == 64)
        assert_equals(18, value['row'][missing])
        assert_equals(64, value['col'][missing])
        assert_equals('N/A', value['id'][missing])
        assert_true(np.isnan(value['ref_signal'][missing]))
        assert_true(np.isnan(value['test_bg'][missing]))
        assert_equals(cghr.INVALID_INT, value['chromosome'][missing])

        # Chromosomes convertions
        Xprobe = (value['id'] == 'A_14_P119856')
        assert_equals(23, value['chromosome'][Xprobe]) # X

    def test_toarray_fields(self):
        value = self.reader.toarray(fields=('row', 'col', 'chromosome'))
        assert_equals(45220, len(value))
        assert_equals(3, len(value[0]))

    def test_toarray_order(self):
        value = self.reader.toarray(order=('chromosome', 'start_base'))
        assert_equals(45220, len(value))
        assert_equals(11, len(value[0]))

        assert_equals(cghr.INVALID_INT, value[0]['chromosome'])
        assert_equals(24, value[-1]['chromosome'])

        # Order with filter
        value = self.reader.toarray(fields=('chromosome',),
                                    order=('chromosome', 'start_base'))
        assert_equals(cghr.INVALID_INT, value[0]['chromosome'])
        assert_equals(24, value[-1]['chromosome'])
        assert_raises(ValueError, value.__getitem__, 'id')

    def test_lazyness(self):
        value1 = self.reader.toarray()
        assert_not_equals(cghr.INVALID_STRING, value1['id'][0])

        value1['id'][0] = cghr.INVALID_STRING
        assert_equals(cghr.INVALID_STRING, value1['id'][0])

        value2 = self.reader.toarray()
        assert_equals(cghr.INVALID_STRING, value2['id'][0])

        value2 = self.reader.toarray(fields=('id',))
        assert_equals(cghr.INVALID_STRING, value2['id'][0])

    def test_swapping(self):
        reader2 = AgilentReader(os.path.join(self.par_dir, 'test_agilent.txt'),
                                test_channel='g')

        value1 = self.reader.toarray()
        value2 = reader2.toarray()

        ref1 = value1['ref_signal'][value1['valid']]
        refbg1 = value1['ref_bg'][value1['valid']]
        test1 = value1['test_signal'][value1['valid']]
        testbg1 = value1['test_bg'][value1['valid']]

        ref2 = value2['ref_signal'][value1['valid']]
        refbg2 = value2['ref_bg'][value1['valid']]
        test2 = value2['test_signal'][value1['valid']]
        testbg2 = value2['test_bg'][value1['valid']]

        assert_true(np.allclose(ref1, test2))
        assert_true(np.allclose(test1, ref2))
        assert_true(np.allclose(refbg1, testbg2))
        assert_true(np.allclose(testbg1, refbg2))
