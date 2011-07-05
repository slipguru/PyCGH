from nose.tools import *
import os

import numpy as np

from pycgh import NimblegenCGH

PAR_DIR = os.path.split(os.path.abspath(__file__))[0]

class TestNimblegenCGH(object):
    def setup(self):
        self.test_path = os.path.join(PAR_DIR, 'test_nimblegen.pair')
        self.ref_path = self.test_path

    def test_access_features(self):
        aCGH = NimblegenCGH.load(self.test_path, self.ref_path)
        assert_equals(24433, len(aCGH['id']))
        assert_equals(24466, len(aCGH.unfiltered('id')))

        assert_raises(ValueError, aCGH.__getitem__, 'fake')

    def test_features_keys(self):
        aCGH = NimblegenCGH.load(self.test_path, self.ref_path)
        assert_equals(9, len(aCGH.names))

    def test_data(self):
        aCGH = NimblegenCGH.load(self.test_path, self.ref_path)

        # Find present data
        present = np.logical_and(aCGH['row'] == 275, aCGH['col'] == 261)
        assert_equals(275, aCGH['row'][present])
        assert_equals(261, aCGH['col'][present])
        assert_equals('CHR0100P000058718', aCGH['id'][present])
        assert_equals(1, aCGH['chromosome'][present])
        assert_equals(58718, aCGH['start_base'][present])
        #assert_equals(11455291, aCGH['end_base'][present])
        assert_true(~np.isnan(aCGH['reference_signal'][present]))

        ## Find missing data
        #missing = np.logical_and(aCGH.unfiltered('row') == 19,
                                 #aCGH.unfiltered('col') == 64)
        #assert_equals(18, aCGH.unfiltered('row')[missing])
        #assert_equals(64, aCGH.unfiltered('col')[missing])
        #assert_equals('N/A', aCGH.unfiltered('id')[missing])
        #assert_true(np.isnan(aCGH.unfiltered('reference_signal')[missing]))
        #assert_equals(NimblegenCGH.INVALID_INT, aCGH.unfiltered('chromosome')[missing])
        #
        # Chromosomes convertions
        Xprobe = (aCGH.unfiltered('id') == 'CHRX00P049585234')
        assert_equals(23, aCGH.unfiltered('chromosome')[Xprobe]) # X
