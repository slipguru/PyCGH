import os

import numpy as np
from numpy.testing import *

from ..datatypes.arraycgh import ArrayCGH
from ..readers import agilent

PAR_PATH = os.path.split(os.path.abspath(__file__))[0]
SAMPLE_PATH = os.path.join(PAR_PATH, 'agilent_sample.txt')

def test_access_features():
    aCGH = agilent(SAMPLE_PATH)
    assert_equal(4460, len(aCGH))
    assert_equal(4460, len(aCGH['id']))
    assert_equal(4233, len(aCGH.F['id']))

    assert_raises(ValueError, aCGH.__getitem__, 'fake')

def test_features_keys():
    aCGH = agilent(SAMPLE_PATH)
    assert_equal(9, len(aCGH.names))

@dec.slow
def test_fill_missing():
    aCGH = agilent(SAMPLE_PATH, fill_missings=True)
    assert_equal(45220, len(aCGH))

def test_qc_cleaning():
    aCGH = agilent(SAMPLE_PATH, qc_masking=False)
    assert_equal(4460, len(aCGH))
    assert_equal(227, aCGH['mask'].sum())

    aCGH = agilent(SAMPLE_PATH, qc_masking=True)
    assert_equal(4460, len(aCGH))
    assert_equal(227 + 3, aCGH['mask'].sum())

@dec.slow
def test_qc_missings():
    aCGH = agilent(SAMPLE_PATH, fill_missings=True, qc_masking=True)
    assert_equal(45220, len(aCGH))
    assert_equal(227 + 3 + (45220 - 4460), aCGH['mask'].sum())

@dec.slow
def test_content():
    aCGH = agilent(SAMPLE_PATH, fill_missings=True)

    # Find present data
    present = np.logical_and(aCGH['row'] == 19, aCGH['col'] == 64)
    assert_equal(19, aCGH['row'][present])
    assert_equal(64, aCGH['col'][present])
    assert_equal('A_14_P102318', aCGH['id'][present])
    assert_equal(19, aCGH['chromosome'][present])
    assert_equal(11455241, aCGH['start_base'][present])
    assert_equal(11455291, aCGH['end_base'][present])
    assert_(~np.isnan(aCGH['reference_signal'][present]))

    # Find missing data
    missing = np.logical_and(aCGH['row'] == 18,
                             aCGH['col'] == 64)
    assert_equal(18, aCGH['row'][missing])
    assert_equal(64, aCGH['col'][missing])
    assert_equal('--', aCGH['id'][missing])
    assert_(np.isnan(aCGH['reference_signal'][missing]))
    assert_equal(ArrayCGH.MISSING_INT, aCGH['chromosome'][missing])
    
    # Chromosomes convertions
    Xprobe = (aCGH['id'] == 'A_14_P119856')
    assert_equal(23, aCGH['chromosome'][Xprobe]) # X

def test_swapping():
    aCGH_1 = agilent(SAMPLE_PATH, test_channel='r')
    aCGH_2 = agilent(SAMPLE_PATH, test_channel='g')

    ref1 = aCGH_1['reference_signal']
    test1 = aCGH_1['test_signal']

    ref2 = aCGH_2['reference_signal']
    test2 = aCGH_2['test_signal']

    assert_equal(ref1, test2)
    assert_equal(test1, ref2)

@dec.slow
def test_masked():
    aCGH_1 = agilent(SAMPLE_PATH, fill_missings=True)
    aCGH_2 = agilent(SAMPLE_PATH, fill_missings=False)

    mask_1 = aCGH_1['mask']
    mask_2 = aCGH_2['mask']
    assert_equal(aCGH_1['reference_signal'][~mask_1],
                 aCGH_2['reference_signal'][~mask_2])
    assert_equal(aCGH_1.F['reference_signal'], aCGH_2.F['reference_signal'])

    # The dimension is different
    assert_(len(aCGH_1.M['reference_signal']) > len(aCGH_2.M['reference_signal']))

    # But not if compressed
    assert_equal(aCGH_1.M['reference_signal'].compressed(),
                 aCGH_2.M['reference_signal'].compressed())

# -- Remapping ----------------------------------------------------------------
from cStringIO import StringIO

RELEASE_FILE = """\
818	chr1	30636423	30636483	A_14_P103951
1935	chr4	177036969	177037029	A_14_P124668
1191	chr8	79514725	79514785	A_14_P102878
"""

def test_release_mapping():
    aCGH = agilent(SAMPLE_PATH)
    index = (aCGH['id'] == 'A_14_P131095')
    assert_equal(6, aCGH['chromosome'][index])
    assert_equal(1259769, aCGH['start_base'][index])
    assert_equal(1259826, aCGH['end_base'][index])

    # Manual remapping of 3 probes
    aCGH = agilent(SAMPLE_PATH, ucsc_mapping = StringIO(RELEASE_FILE))
    assert_equal(3, len(aCGH.F['id']))
    
    index = (aCGH['id'] == 'A_14_P103951')    
    assert_equal(1, aCGH['chromosome'][index])
    assert_equal(30636424, aCGH['start_base'][index])  # UCSC format shift
    assert_equal(30636483, aCGH['end_base'][index])
