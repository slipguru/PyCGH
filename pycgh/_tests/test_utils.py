import os
from nose.tools import *

from pycgh.utils import CytoBands, probes_average, LabeledMatrix

PAR_DIR = os.path.split(os.path.abspath(__file__))[0]



class TestProbesAverage(object):
    def setup(self):
        self.ids =    [1,   2,  3,  1,  1,  4,  5]
        self.values = [10, 10, 10, 20, 300, 10, 10]

    def test_average(self):
        import numpy as np

        out = probes_average(self.ids, self.values, np.mean)
        assert_equals(110, out[1])

        out = probes_average(self.ids, self.values, np.median)
        assert_equals(20, out[1])

        # Unique values
        assert_equals(10, out[2])
        assert_equals(10, out[3])
        assert_equals(10, out[4])
        assert_equals(10, out[5])

