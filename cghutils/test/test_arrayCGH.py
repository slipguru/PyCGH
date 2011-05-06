import os

from nose.tools import *

import itertools as it
import numpy as np
from cghutils import ArrayCGH

PAR_DIR = os.path.split(os.path.abspath(__file__))[0]
ROW_NUM = COL_NUM = 10

class TestArrayCGH(object):

    def setup(self):
        row, col =  zip(*it.product(range(1, ROW_NUM+1),
                                    range(1, COL_NUM+1)))
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
        assert_equals(9, len(aCGH.names))

    def test_missing_mandatory(self):
        assert_raises(ValueError, ArrayCGH, self.input[:-1])

    def test_optional_inputs(self):
        mask = [False]*(ROW_NUM*COL_NUM)
        aCGH = ArrayCGH(self.input, mask=mask)
        assert_equals(100, len(aCGH))
        assert_equals(9, len(aCGH.names))

        optional = [30]*(ROW_NUM*COL_NUM)
        aCGH = ArrayCGH(self.input, optional_signal=optional)
        assert_equals(10, len(aCGH.names))

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

    def test_masked(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(self.input, mask=mask)

        assert_equals(ROW_NUM*COL_NUM, len(aCGH['id']))
        assert_equals((ROW_NUM*COL_NUM) - 10, len(aCGH.filtered('id')))

        assert_equals((ROW_NUM*COL_NUM), len(aCGH.masked('id')))
        assert_equals((ROW_NUM*COL_NUM) - 10, len(aCGH.masked('id').compressed()))

    def test_adding(self):
        aCGH = ArrayCGH(self.input)

        ratio = np.log2(aCGH['test_signal']/aCGH['reference_signal'])
        assert_equals(ROW_NUM*COL_NUM, len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.allclose(ratio, aCGH['ratio']))
        assert_true('ratio' in aCGH.names)

    def test_adding_filtered(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(self.input, mask=mask)

        ratio = np.log2(aCGH.filtered('test_signal')/aCGH.filtered('reference_signal'))
        assert_equals((ROW_NUM*COL_NUM) - 10, len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.allclose(ratio, aCGH.filtered('ratio')))
        assert_true('ratio' in aCGH.names)

        assert_raises(ValueError, aCGH.__setitem__, 'foo', ratio[:-5]) #shorter
        assert_raises(ValueError, aCGH.__setitem__, 'foo', np.r_[ratio, ratio]) #longer

    def test_adding_masked(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(self.input, mask=mask)

        ratio = np.log2(aCGH.masked('test_signal')/aCGH.masked('reference_signal'))
        assert_equals((ROW_NUM*COL_NUM), len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.ma.allclose(ratio, aCGH.masked('ratio'))) # in "ma" sense
        assert_true('ratio' in aCGH.names)

    def test_update(self):
        aCGH = ArrayCGH(self.input)
        aCGH['mask'] = aCGH['mask']
        assert_true(False)
        #TODO: TEST TEST TEST


class TestArrayCGHIO(object):

    def test_loading_file(self):
        aCGH = ArrayCGH.load(os.path.join(PAR_DIR, 'test_acgh.txt'))
        assert_true(np.all(np.array([1, 1, 2, 2]) == aCGH['row']))

    def test_loading_fields(self):
        aCGH = ArrayCGH.load(os.path.join(PAR_DIR, 'test_acgh2.txt'),
                             fields={'id':'ProbeID',
                                     'reference_signal':'Ref',
                                     'test_signal': 'Test',
                                     'flag': 'myflag'}) # Name change
        assert_true(np.all(np.array([0, 0, 0, 1]) == aCGH['mask']))

    def test_save(self):
        import tempfile
        out = os.path.join(tempfile.gettempdir(), 'acgh.txt')

        aCGH = ArrayCGH.load(os.path.join(PAR_DIR, 'test_acgh.txt'),
                             fields={'orig_mask':'mask', #renaming
                                     'mask':'myflag'})
        aCGH.save(out)
        aCGH2 = ArrayCGH.load(out)

        assert_equals(len(aCGH), len(aCGH2))
        assert_equals(aCGH.names, aCGH2.names)
        for a, b in zip(aCGH, aCGH):
            assert_equals(a[0], b[0])
