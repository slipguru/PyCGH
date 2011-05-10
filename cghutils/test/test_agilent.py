from nose.tools import *
import os

import numpy as np

from cghutils import AgilentCGH

PAR_DIR = os.path.split(os.path.abspath(__file__))[0]

class TestAgilentCGH(object):
    def setup(self):
        self.path = os.path.join(PAR_DIR, 'test_agilent.txt')

    def test_access_features(self):
        aCGH = AgilentCGH.load(self.path)
        assert_equals(4460, len(aCGH))
        assert_equals(4460, len(aCGH.unfiltered('id')))

        assert_raises(ValueError, aCGH.__getitem__, 'fake')

    def test_features_keys(self):
        aCGH = AgilentCGH.load(self.path)
        assert_equals(9, len(aCGH.names))

    def test_fill_missing(self):
        aCGH = AgilentCGH.load(self.path, fill_missings=True)
        assert_equals(45220, len(aCGH))

    def test_qc_cleaning(self):
        aCGH = AgilentCGH.load(self.path, qc_masking=False)
        assert_equals(4460, len(aCGH))
        assert_equals(227, aCGH.unfiltered('mask').sum())

        aCGH = AgilentCGH.load(self.path, qc_masking=True)
        assert_equals(4460, len(aCGH))
        assert_equals(227+3, aCGH.unfiltered('mask').sum())

    def test_qc_missings(self):
        aCGH = AgilentCGH.load(self.path, fill_missings=True, qc_masking=True)
        assert_equals(45220, len(aCGH))
        assert_equals(227+3 + (45220-4460), aCGH.unfiltered('mask').sum())

    def test_data(self):
        aCGH = AgilentCGH.load(self.path, fill_missings=True)

        # Find present data
        present = np.logical_and(aCGH['row'] == 19, aCGH['col'] == 64)
        assert_equals(19, aCGH['row'][present])
        assert_equals(64, aCGH['col'][present])
        assert_equals('A_14_P102318', aCGH['id'][present])
        assert_equals(19, aCGH['chromosome'][present])
        assert_equals(11455241, aCGH['start_base'][present])
        assert_equals(11455291, aCGH['end_base'][present])
        assert_true(~np.isnan(aCGH['reference_signal'][present]))

        # Find missing data
        missing = np.logical_and(aCGH.unfiltered('row') == 18,
                                 aCGH.unfiltered('col') == 64)
        assert_equals(18, aCGH.unfiltered('row')[missing])
        assert_equals(64, aCGH.unfiltered('col')[missing])
        assert_equals('N/A', aCGH.unfiltered('id')[missing])
        assert_true(np.isnan(aCGH.unfiltered('reference_signal')[missing]))
        assert_equals(AgilentCGH.INVALID_INT, aCGH.unfiltered('chromosome')[missing])

        # Chromosomes convertions
        Xprobe = (aCGH.unfiltered('id') == 'A_14_P119856')
        assert_equals(23, aCGH.unfiltered('chromosome')[Xprobe]) # X

    def test_swapping(self):
        aCGH_1 = AgilentCGH.load(self.path, test_channel='r')
        aCGH_2 = AgilentCGH.load(self.path, test_channel='g')

        ref1 = aCGH_1['reference_signal']
        test1 = aCGH_1['test_signal']

        ref2 = aCGH_2['reference_signal']
        test2 = aCGH_2['test_signal']

        assert_true(np.allclose(ref1, test2))
        assert_true(np.allclose(test1, ref2))

    def test_masked(self):
        aCGH_1 = AgilentCGH.load(self.path, fill_missings=True)
        aCGH_2 = AgilentCGH.load(self.path, fill_missings=False)

        mask_1 = aCGH_1['mask']
        mask_2 = aCGH_2['mask']
        assert_true(np.allclose(aCGH_1.unfiltered('reference_signal')[~mask_1],
                                aCGH_2.unfiltered('reference_signal')[~mask_2]))

        assert_true(np.allclose(aCGH_1['reference_signal'],
                                aCGH_2['reference_signal']))

        # The dimension is different
        assert_false(np.allclose(aCGH_1.masked('reference_signal'),
                                 aCGH_2.masked('reference_signal')))

        # But true if compressed
        assert_true(np.allclose(aCGH_1.masked('reference_signal').compressed(),
                                 aCGH_2.masked('reference_signal').compressed()))
