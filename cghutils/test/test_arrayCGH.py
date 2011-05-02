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

        self.input = [['Probe%d' % x for x in range(1, len(row)+1)],    #ID
                      list(row),                                        #Row
                      list(col),                                        #Col
                      [10]*len(row),                                    #Ref
                      [20]*len(row),                                    #Test
                      (range(1, 25)*5)[:len(row)],                      #Chr
                      range(1, len(row)+1),                             #Start
                      [x+10 for x in range(1, len(row)+1)]]             #End

    def test_simple_creation(self):
        aCGh = ArrayCGH(self.input)

    def test_optional_inputs(self):
        mask = [True]*(self.ROW_NUM*self.COL_NUM)
        aCGH = ArrayCGH(self.input, mask=mask)

        optional = [30]*(self.ROW_NUM*self.COL_NUM)
        aCGH = ArrayCGH(self.input, optional_signal=optional)

        aCGH = ArrayCGH(self.input, mask=mask, optional_signal=optional)

        # name duplication
        assert_raises(ValueError, ArrayCGH, self.input, chromosome=optional)
