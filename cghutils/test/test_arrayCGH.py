import os

from nose.tools import *

import itertools as it
import numpy as np
from cghutils import ArrayCGH

class TestArrayCGH(object):

    def setup(self):
        self.PAR_DIR = os.path.split(os.path.abspath(__file__))[0]
        self.ROW_NUM = self.COL_NUM = 10
        row, col =  zip(*it.product(range(1, self.ROW_NUM+1),
                                    range(1, self.COL_NUM+1)))
        self.input = [['Probe%d' % x for x in range(1, len(row)+1)],    #ID
                      list(row),                                        #Row
                      list(col),                                        #Col
                      [10]*len(row),                                    #Ref
                      [20]*len(row),                                    #Test
                      (range(1, 25)*5)[:len(row)],                      #Chr
                      range(1, len(row)+1),                             #Start
                      [x+10 for x in range(1, len(row)+1)]]             #End

    def test_simple_creation(self):
        aCGH = ArrayCGH(self.input)
        assert_equals(100, len(aCGH))
        assert_equals(8, len(aCGH.names))

    def test_optional_inputs(self):
        mask = [True]*(self.ROW_NUM*self.COL_NUM)
        aCGH = ArrayCGH(self.input, mask=mask)
        assert_equals(100, len(aCGH))
        assert_equals(9, len(aCGH.names))

        optional = [30]*(self.ROW_NUM*self.COL_NUM)
        aCGH = ArrayCGH(self.input, optional_signal=optional)
        assert_equals(9, len(aCGH.names))

        aCGH = ArrayCGH(self.input, mask=mask, optional_signal=optional)
        assert_equals(10, len(aCGH.names))

        # name duplication
        assert_raises(ValueError, ArrayCGH, self.input, chromosome=optional)

    def test_reading(self):
        aCGH = ArrayCGH(self.input)

        id = aCGH['id']
        assert_true(np.equal(self.input[0], id))

        submatrix = aCGH[['id', 'chromosome']]
        assert_true(np.all(self.input[0] == submatrix['id']))
        assert_true(np.all(self.input[5] == submatrix['chromosome']))

    def test_loading_file(self):
        aCGH = ArrayCGH.load(os.path.join(self.PAR_DIR, 'test_acgh.txt'))
        assert_true(np.all(np.array([1, 1, 2, 2]) == aCGH['row']))

    def test_loading_fields(self):
        aCGH = ArrayCGH.load(os.path.join(self.PAR_DIR, 'test_acgh2.txt'),
                             fields={'id':'ProbeID',
                                     'reference_signal':'Ref',
                                     'test_signal': 'Test',
                                     'flag': 'myflag'}) # Name change
        assert_true(np.all(np.array([0, 0, 0, 1]) == aCGH['mask']))

    def test_save(self):
        import tempfile
        out = os.path.join(tempfile.gettempdir(), 'acgh.txt')

        aCGH = ArrayCGH.load(os.path.join(self.PAR_DIR, 'test_acgh.txt'),
                             fields={'orig_mask':'mask', #renaming
                                     'mask':'myflag'})
        aCGH.save(out)
        aCGH2 = ArrayCGH.load(out)

        assert_equals(len(aCGH), len(aCGH2))
        assert_equals(aCGH.names, aCGH2.names)
        for a, b in zip(aCGH, aCGH):
            assert_equals(a[0], b[0])
