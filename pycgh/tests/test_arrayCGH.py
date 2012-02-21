import itertools as it
from cStringIO import StringIO

import numpy as np
from numpy.testing.utils import *

from ..datatypes.arraycgh import ArrayCGH

def setup(module):
    ROW_NUM = COL_NUM = 10
    row, col =  zip(*it.product(range(1, ROW_NUM+1),
                                range(1, COL_NUM+1)))

    # 10x10 aCGH support with 100 probes
    input = (['Probe%d' % x for x in range(1, len(row)+1)],         #ID
             list(row),                                             #Row
             list(col),                                             #Col
             [10] * (ROW_NUM * COL_NUM),                            #Ref
             [20] * (ROW_NUM * COL_NUM),                            #Test
             (range(1, 25) * 5)[:(ROW_NUM * COL_NUM)],              #Chr
             range(1, (ROW_NUM * COL_NUM) +1 ),                     #Start
             [x + 10 for x in range(1, (ROW_NUM * COL_NUM) + 1)])   #End

    # Setup global variables
    module.input = input
    module.ROW_NUM = ROW_NUM
    module.COL_NUM = COL_NUM

def test_simple_creation():
    aCGH = ArrayCGH(*input)
    assert_equal(100, len(aCGH))
    assert_equal(9, len(aCGH.names))
    
def test_missing_mandatory():
    assert_raises(TypeError, ArrayCGH, input[:-1])

def test_optional_inputs():
    mask = [False] * (ROW_NUM * COL_NUM)
    aCGH = ArrayCGH(*input, mask=mask) # as default
    assert_equal(100, len(aCGH))
    assert_equal(9, len(aCGH.names))

    optional = [30] * (ROW_NUM * COL_NUM)
    aCGH = ArrayCGH(*input, optional_signal=optional)
    assert_equal(10, len(aCGH.names))

    # Mask and other optional arguments
    aCGH = ArrayCGH(*input, mask=mask, optional_signal=optional)
    assert_equal(10, len(aCGH.names))

    # name duplication
    assert_raises(TypeError, ArrayCGH, input, chromosome=optional)

def test_reading():
    aCGH = ArrayCGH(*input)

    id = aCGH['id']
    assert_equal(input[0], id)

    submatrix = aCGH[['id', 'chromosome']]
    assert_equal(input[0], submatrix['id'])
    assert_equal(input[5], submatrix['chromosome'])

def test_masked():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 masked values
    aCGH = ArrayCGH(*input, mask=mask)

    assert_equal(ROW_NUM * COL_NUM, len(aCGH.unfiltered('id')))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH['id']))

    assert_equal((ROW_NUM * COL_NUM), len(aCGH.masked('id')))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.masked('id').compressed()))


class _TestArrayCGH(object):    



    def test_adding_unfiltered(self):
        aCGH = ArrayCGH(*self.input)

        ratio = np.log2(aCGH.unfiltered('test_signal')/aCGH.unfiltered('reference_signal'))
        assert_equals(ROW_NUM*COL_NUM, len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.allclose(ratio, aCGH.unfiltered('ratio')))
        assert_true('ratio' in aCGH.names)

    def test_adding(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(*self.input, mask=mask)

        ratio = np.log2(aCGH['test_signal']/aCGH['reference_signal'])
        assert_equals((ROW_NUM*COL_NUM) - 10, len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.allclose(ratio, aCGH['ratio']))
        assert_true('ratio' in aCGH.names)

        assert_raises(ValueError, aCGH.__setitem__, 'foo', ratio[:-5]) #shorter
        assert_raises(ValueError, aCGH.__setitem__, 'foo', np.r_[ratio, ratio]) #longer

    def test_adding_masked(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(*self.input, mask=mask)

        ratio = np.log2(aCGH.masked('test_signal')/aCGH.masked('reference_signal'))
        assert_equals((ROW_NUM*COL_NUM), len(ratio))

        aCGH['ratio'] = ratio
        assert_true(np.ma.allclose(ratio, aCGH.masked('ratio'))) # in "ma" sense
        assert_true('ratio' in aCGH.names)

    def test_update_unfiltered(self):
        aCGH = ArrayCGH(*self.input)
        new_mask = np.array([True]*(ROW_NUM*COL_NUM))
        aCGH['mask'] = new_mask

        assert_true(np.allclose(new_mask, aCGH.unfiltered('mask')))
        assert_equals(0, len(aCGH['id']))

    def test_update(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(*self.input, mask=mask)

        filt_test = aCGH['test_signal']
        filt_test += 1.
        assert_equals((~mask).sum(), len(filt_test))

        assert_false(np.allclose(filt_test, aCGH['test_signal']))

        aCGH['test_signal'] = filt_test
        assert_true(np.allclose(filt_test, aCGH['test_signal']))

        assert_raises(ValueError, aCGH.__setitem__, 'test_signal', filt_test[:-5]) #shorter
        assert_raises(ValueError, aCGH.__setitem__, 'test_signal', np.r_[filt_test, filt_test]) #longer

    def test_update_masked(self):
        mask = np.array([False]*(ROW_NUM*COL_NUM))
        mask[0:10] = True # 10 values masked
        aCGH = ArrayCGH(*self.input, mask=mask)

        mask_test = aCGH.masked('test_signal', copy=True) # Force copying
        mask_test += 1.
        assert_equals((ROW_NUM*COL_NUM), len(mask_test))
        assert_equals((~mask).sum(), len(mask_test.compressed()))
        assert_false(np.ma.allclose(mask_test, aCGH.masked('test_signal')))

        aCGH['test_signal'] = mask_test
        assert_true(np.ma.allclose(mask_test, aCGH.masked('test_signal')))

    def test_order(self):
        # Call this in-place ordering may be useful before saving the data
        # or before the estraction
        aCGH = ArrayCGH(*self.input)

        chr = aCGH['chromosome'].copy()
        idx = np.argsort(chr)

        aCGH.sort('chromosome')
        assert_true(np.allclose(chr[idx], aCGH['chromosome']))




aCGHContent = """\
ProbeId, col, row, Ref, Test, chromosome, start_base, end_base, mask, myflag
Probe01, 1, 1, 10, 20, 1, 10, 20, 0, 1
Probe02, 2, 1, 10, 20, 2, 10, 20, 0, 1
Probe 03, 1, 2, 10, 20, 3, 10, 20, 0, 1
Control , 2, 2, 10, 20, nan, nan, nan, 1, 1
"""

"""\
id, col, row, Reference_signal, test_signal, chromosome, start_base, end_base, mask, myflag
Probe01, 1, 1, 10, 10, 1, 10, 10.123456, 0, 1
Probe02, 2, 1, 10, 10, 2, 10, 20, 0, 1
Probe 03, 1, 2, 10, 10, 3, 10000.1039478, 20, 0, 1
Control, 2, 2, 10, 10, nan, nan, nan, 1, 1
"""


class _TestArrayCGHIO(object):

    def test_loading_file(self):
        aCGH = ArrayCGH.load(os.path.join(PAR_DIR, 'test_acgh.txt'))
        assert_true(np.all(np.array([1, 1, 2, 2]) == aCGH.unfiltered('row')))

    def test_loading_fields(self):
        aCGH = ArrayCGH.load(os.path.join(PAR_DIR, 'test_acgh2.txt'),
                             fields={'id':'ProbeID',
                                     'reference_signal':'Ref',
                                     'test_signal': 'Test',
                                     'flag': 'myflag'}) # Name change
        assert_true(np.all(np.array([0, 0, 0, 1]) == aCGH.unfiltered('mask')))

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
