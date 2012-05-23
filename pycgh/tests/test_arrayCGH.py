import itertools as it
from cStringIO import StringIO

import numpy as np
from numpy.testing.utils import *

from ..datatypes.arraycgh import ArrayCGH

# Setup global variables ------------------------------------------------------
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
# -----------------------------------------------------------------------------

def test_types():
    aCGH = ArrayCGH(id=['a', 'b'], row=[1., 1], col=[1., 2.],
                    reference_signal=[1, '2'], test_signal=[1., 2],
                    chromosome=[1, 2], start_base=[1, 1], end_base=[2, 2])

    assert_equal(['a', 'b'], aCGH['id'])
    assert_equal([1, 1], aCGH['row'])
    assert_equal([1, 2], aCGH['col'])
    assert_equal([1., 2.], aCGH['reference_signal'])
    assert_equal([1., 2.], aCGH['test_signal'])
    assert_equal([1., 2.], aCGH['chromosome'])
    assert_equal([1, 1], aCGH['start_base'])
    assert_equal([2, 2], aCGH['end_base'])

    assert_equal(np.dtype('|S1'), aCGH['id'].dtype)
    assert_equal(int, aCGH['row'].dtype)
    assert_equal(int, aCGH['col'].dtype)
    assert_equal(float, aCGH['reference_signal'].dtype)
    assert_equal(float, aCGH['test_signal'].dtype)
    assert_equal(int, aCGH['chromosome'].dtype)
    assert_equal(int, aCGH['start_base'].dtype)
    assert_equal(int, aCGH['end_base'].dtype)

def test_chromosome_conversion():
    aCGH = ArrayCGH(id=['a', 'b'], row=[1., 1], col=[1., 2.],
                    reference_signal=[1, '2'], test_signal=[1., 2],
                    chromosome=['X', 'Y'], start_base=[1, 1], end_base=[2, 2])

    assert_equal([23, 24], aCGH['chromosome'])
    assert_equal(int, aCGH['chromosome'].dtype)

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

    pairs = sorted((a, b) for a, b in zip(input[5], input[6]))
    chr, sb = zip(*pairs)

    submatrix = aCGH[['chromosome', 'start_base']]
    assert_equal(chr, submatrix['chromosome'])
    assert_equal(sb, submatrix['start_base'])

def test_masked_filtered():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 masked values
    aCGH = ArrayCGH(*input, mask=mask)

    # Single col
    assert_equal(ROW_NUM * COL_NUM, len(aCGH['id']))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.filtered('id')))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.F['id'])) # alternative

    # Multiple col
    assert_equal(ROW_NUM * COL_NUM, len(aCGH[['id', 'mask']]))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.filtered(['id', 'mask'])))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.F[['id', 'mask']]))

    # Masked and different number of columns
    assert_equal((ROW_NUM * COL_NUM), len(aCGH.masked(['id', 'mask'])))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.masked('id').compressed()))

    assert_equal((ROW_NUM * COL_NUM), len(aCGH.M[['id', 'mask']])) # alt
    assert_equal((ROW_NUM * COL_NUM) - 10, len(aCGH.M['id'].compressed()))

def test_adding():
    aCGH = ArrayCGH(*input)

    ratio = np.log2(aCGH['test_signal']/aCGH['reference_signal'])
    assert_equal(ROW_NUM * COL_NUM, len(ratio))

    aCGH['ratio'] = ratio
    assert_equal(ratio, aCGH['ratio'])
    assert_equal(np.compress(~aCGH['mask'], ratio), aCGH.filtered('ratio'))
    assert_equal(np.ma.masked_array(aCGH['ratio'], aCGH['mask']),
                 aCGH.masked('ratio'))
    assert_('ratio' in aCGH.names)

def test_adding_filtered():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 values masked
    aCGH = ArrayCGH(*input, mask=mask)

    ratio = np.log2(aCGH.filtered('test_signal') /
                    aCGH.filtered('reference_signal'))
    assert_equal((ROW_NUM * COL_NUM) - 10, len(ratio))

    aCGH['ratio'] = ratio # automatic filtering
    assert_equal(ratio, aCGH.filtered('ratio'))
    assert_('ratio' in aCGH.names)

    assert_raises(ValueError, aCGH.__setitem__, 'foo', ratio[:-5]) #shorter
    assert_raises(ValueError, aCGH.__setitem__, 'foo', np.r_[ratio, ratio]) #longer

def test_adding_masked():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 values masked
    aCGH = ArrayCGH(*input, mask=mask)

    ratio = np.log2(aCGH.M['test_signal'] /
                    aCGH.M['reference_signal'])
    assert_equal((ROW_NUM * COL_NUM), len(ratio))

    aCGH['ratio'] = ratio # automatic masking
    assert_equal(ratio, aCGH.M['ratio']) # in "ma" sense
    assert_('ratio' in aCGH.names)

def test_update():
    aCGH = ArrayCGH(*input)
    new_mask = np.array([True] * (ROW_NUM * COL_NUM))
    aCGH['mask'] = new_mask

    assert_equal(new_mask, aCGH['mask'])
    assert_equal(0, len(aCGH.filtered('id')))

def test_update():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 values masked
    aCGH = ArrayCGH(*input, mask=mask)

    filt_test = aCGH.F['test_signal'] + 1.
    assert_(not np.allclose(filt_test, aCGH.F['test_signal']))

    aCGH['test_signal'] = filt_test
    assert_equal(filt_test, aCGH.F['test_signal'])

    assert_raises(ValueError, aCGH.__setitem__, 'test_signal',
                  filt_test[:-5]) #shorter
    assert_raises(ValueError, aCGH.__setitem__, 'test_signal',
                  np.r_[filt_test, filt_test]) #longer

def test_update_masked():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 values masked
    aCGH = ArrayCGH(*input, mask=mask)

    mask_test = aCGH.masked('test_signal', copy=True) + 1# Force copying
    assert_equal((ROW_NUM * COL_NUM), len(mask_test))
    assert_(not np.ma.allclose(mask_test, aCGH.M['test_signal']))

    aCGH['test_signal'] = mask_test
    assert_equal(mask_test, aCGH.M['test_signal'])

def test_default_order():
    # Standard input is not sorted by chromosome and bases
    aCGH = ArrayCGH(*input)
    input_pairs = [tuple(x) for x in aCGH[['chromosome', 'start_base']]]

    input_pairs.sort()
    assert_equal(input_pairs,
                 [tuple(x) for x in aCGH[['chromosome', 'start_base']]])

def test_shrinking():
    mask = np.array([False] * (ROW_NUM * COL_NUM))
    mask[0:10] = True # 10 values masked
    aCGH = ArrayCGH(*input, mask=mask)

    assert_equal(100, len(aCGH))
    assert_equal(90, aCGH.size)

    aCGH.shrink()

    assert_equal(aCGH.size, len(aCGH))
    assert_equal(90, aCGH.size)
    assert_equal(aCGH['id'], aCGH.F['id'])
    assert_equal(aCGH['id'], aCGH.M['id'])


## Testing IO -----------------------------------------------------------------

aCGHContent = """\
id, col, row, Reference_signal, test_signal, chromosome, start_base, end_base, mask, myflag
Probe01, 1, 1, 10., 10, 1, 10, 10.123456, 0, 1
Probe02, 2, 1, 10., 10, 2, 10, 20, 0, 1
Probe 03, 1, 2, 10., 10, X, 10000.1039478, 20, 0, 1
######
# This is a Comment
######
Control, 2, 2, --, 10, --, --, --, 1, 1
"""

aCGHContentNotStandard = """\
ProbeId, col, row, Ref, Test, chromosome, start_base, end_base, mask, myflag
Probe01, 1, 1, 10, 20, 1, 10, 20, 0, 1
Probe02, 2, 1, 10, 20, 2, 10, 20, 0, 1
######
# This is a Comment
######
Probe 03, 1, 2, 10, 20, 3, 10, 20, 0, 1
Control , 2, 2, 10, 20, --, --, --, 1, 1
# This is a Comment
"""

def test_loading_file():
    aCGH = ArrayCGH.loadtxt(StringIO(aCGHContent))
    assert_equal([2, 1, 1, 2], aCGH['row'])
    assert_equal(['Control', 'Probe01', 'Probe02', 'Probe 03'], aCGH['id'])

    assert_equal([-1, 1, 2, 23], aCGH['chromosome'])
    assert_equal([np.nan, 10., 10., 10.], aCGH['reference_signal'])

    # Automatic truncated
    assert_equal([-1, 10, 10, 10000], aCGH['start_base'])
    assert_equal([-1, 10, 20, 20], aCGH['end_base'])

def test_loading_fields():
    aCGH = ArrayCGH.loadtxt(StringIO(aCGHContentNotStandard),
                            fields={'id':'ProbeID',
                                    'reference_signal':'Ref',
                                    'test_signal': 'Test',
                                    'flag': 'myflag'}) # Renaming
    assert_equal([1, 0, 0, 0], aCGH['mask'])
    assert_equal(['Control', 'Probe01', 'Probe02', 'Probe 03'], aCGH['id'])

def test_savetxt():
    aCGH = ArrayCGH.loadtxt(StringIO(aCGHContent),
                         fields={'orig_mask':'mask', #renaming
                                 'mask':'myflag'})

    # Saving
    out = StringIO()
    aCGH.savetxt(out)

    # Reading saved file
    aCGH2 = ArrayCGH.loadtxt(StringIO(out.getvalue()))

    assert_equal(len(aCGH), len(aCGH2))
    assert_equal(aCGH.names, aCGH2.names)
    assert_('orig_mask' in aCGH2.names)
    assert_(not 'myflag' in aCGH2.names)

    for n in aCGH.names:
        assert_equal(aCGH[n], aCGH2[n])

def test_save():
    aCGH = ArrayCGH.loadtxt(StringIO(aCGHContent),
                         fields={'orig_mask':'mask', #renaming
                                 'mask':'myflag'})

    # Saving
    from tempfile import TemporaryFile
    out = TemporaryFile()
    aCGH.save(out)

    # Reading saved file
    out.seek(0)
    aCGH2 = ArrayCGH.load(out)

    assert_equal(len(aCGH), len(aCGH2))
    assert_equal(aCGH.names, aCGH2.names)
    assert_('orig_mask' in aCGH2.names)
    assert_(not 'myflag' in aCGH2.names)

    for n in aCGH.names:
        assert_equal(aCGH[n], aCGH2[n])
