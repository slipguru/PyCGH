from nose.tools import *

import itertools as it
import numpy as np
from cghutils import ArrayCGH

class TestArrayCGH(object):

    def setup(self):
        self.ROW_NUM = self.COL_NUM = 10
        row, col =  zip(*it.product(range(1, self.ROW_NUM+1),
                                    range(1, self.COL_NUM+1)))
        self.input =  {'row': list(row),
                       'col': list(col),
                       'id': ['Probe%d' % x for x in range(1, len(row)+1)],
                       'ref_signal': [10]*len(row),
                       'sample_signal': [20]*len(row),
                       'ref_bg': [1]*len(row),
                       'sample_bg': [2]*len(row),
                       'chromosome': (range(1, 25)*5)[:len(row)],
                       'start_base': range(1, len(row)+1),
                       'end_base': [x+10 for x in range(1, len(row)+1)]
                    }

    def test_subscription(self):
        acgh = ArrayCGH(self.input)
        for k in self.input:
            assert_equals(len(self.input[k]), len(acgh[k]))

    def test_ratios(self):
        acgh = ArrayCGH(self.input)
        ratios = acgh.compute_ratios()

        assert_equals(self.ROW_NUM*self.COL_NUM, len(ratios))

        expected_ratios = np.log2(acgh['sample_signal'] /
                                  acgh['ref_signal'])
        assert_true(np.allclose(expected_ratios, ratios))

    def test_ratios_bg(self):
        acgh = ArrayCGH(self.input)
        ratios = acgh.compute_ratios(bg_subtracted=True)

        expected_ratios = np.log2(
                            (acgh['sample_signal'] - acgh['sample_bg']) /
                            (acgh['ref_signal'] - acgh['ref_bg']))
        assert_true(np.allclose(expected_ratios, ratios))

    def test_ids(self):
        acgh = ArrayCGH(self.input)
        assert_equals(np.dtype('|S8'), acgh['id'].dtype)

    def test_wrong_shape(self):
        for k in self.input:
            input = dict(self.input) # copy
            input[k] = input[k][1:] # removing one element

            try:
                ArrayCGH(self.input)
            except ValueError, e:
                assert_equals("array-shape mismatch in '%s'" % k, str(e))

    def test_filtered_data(self):
        filter = np.array([True]*len(self.input['row']))
        filter[10] = False #masking one element

        acgh = ArrayCGH(self.input, filter=filter)
        for k in self.input:
            assert_equals(len(self.input[k])-1, len(acgh[k].compressed()))

    #def test_ratios_ordered(self):
        #acgh = ArrayCGH(self.input)
        #ratios = acgh.compute_ratios(ordered=True)
