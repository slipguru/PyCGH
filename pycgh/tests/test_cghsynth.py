#import itertools as it
from cStringIO import StringIO

import numpy as np
from numpy.testing.utils import *

from ..utils import ArrayCGHSynth
from ..datatypes.arraycgh import ArrayCGH
from ..datatypes.cytobands import CytoStructure

CHIP_DESIGN = {
    'P01': (1, 1, 200),
    'P02': (1, 300, 400),
    'P03': (1, 500, 600),
    'P04': (1, 700, 800),
    'P05': (2, 1, 200),
}

# UCSC format
CytoFileContent = """\
chr1	0	200	p12	    gpos50
chr1	200	400	p11.2	gneg
chr1	400	800	p11.1	acen
chr2	0	150	q11.2	gneg
"""

NROW = NCOL = 10
NPROBES = NROW * NCOL

def test_default():
    cgh_src = ArrayCGHSynth(geometry=(NROW, NCOL),
                            design=CHIP_DESIGN) # replicas??

    acgh = cgh_src.draw()
    assert_equal(ArrayCGH, type(acgh))
    assert_equal(NPROBES, len(acgh))
    assert_equal(len(CHIP_DESIGN), len(acgh.F['id']))

def test_order():
    cgh_src = ArrayCGHSynth(geometry=(NROW, NCOL),
                            design=CHIP_DESIGN)
    acgh = cgh_src.draw()

    # Sorting by chromosome order
    acgh.sort()
    CHIP_SORTED = sorted(CHIP_DESIGN.items())

    assert_equal([x[0] for x in CHIP_SORTED], acgh.F['id'])
    assert_equal([x[1][0] for x in CHIP_SORTED], acgh.F['chromosome'])
    assert_equal([x[1][1] for x in CHIP_SORTED], acgh.F['start_base'])
    assert_equal([x[1][2] for x in CHIP_SORTED], acgh.F['end_base'])

def test_missing_cytostructure():
    # Non-Empty alteration without cytostructure
    assert_raises(ValueError, ArrayCGHSynth, (NROW, NCOL), CHIP_DESIGN,
                  {'1':[(2, 1.0)]})

def test_parameters_geometry():
    cgh_src = ArrayCGHSynth(design=CHIP_DESIGN, geometry=(NROW, NCOL))
    assert_equal((NROW, NCOL), cgh_src.geometry)

    # Not tuple
    assert_raises(TypeError, ArrayCGHSynth, NROW, CHIP_DESIGN)

    # Bad tuple
    assert_raises(ValueError, ArrayCGHSynth, (-5, 2), CHIP_DESIGN)

    # Truncated
    cgh_src = ArrayCGHSynth(design=CHIP_DESIGN, geometry=(1.2, 7.4))
    assert_equal((1, 7), cgh_src.geometry)

def test_parameters_noise():
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, noise=(0.2, 0.5))
    assert_equal((0.2, 0.5), cgh_src.noise)

    # Default
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN)
    assert_equal((0.1, 0.2), cgh_src.noise)

    # Not tuple
    assert_raises(TypeError, ArrayCGHSynth, (NROW, NCOL), CHIP_DESIGN,
                                            noise=0.2)

    # Sorting
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, noise=(0.5, 0.2))
    assert_equal((0.2, 0.5), cgh_src.noise)

    # Bad tuple
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, noise=(-0.5, 0.2))

def test_parameter_tissue_proportion():
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, tissue_prop=(0.2, 0.4))
    assert_equal((0.2, 0.4), cgh_src.tissue_prop)

    # Default
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN)
    assert_equal((0.3, 0.7), cgh_src.tissue_prop)

    # Not tuple
    assert_raises(TypeError, ArrayCGHSynth, (NROW, NCOL), CHIP_DESIGN,
                                            tissue_prop=0.2)

    # Sorting
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, tissue_prop=(0.5, 0.2))
    assert_equal((0.2, 0.5), cgh_src.tissue_prop)

    # Bad tuple
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, tissue_prop=(-0.5, 0.2))
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, tissue_prop=(0.5, 1.01))
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, tissue_prop=(1.1, 1.))

def test_parameters_dyes():
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, dyes=(50, 100))
    assert_equal((50, 100), cgh_src.dyes)

    # Default
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN)
    assert_equal((100, 500), cgh_src.dyes)

    # Not tuple
    assert_raises(TypeError, ArrayCGHSynth, (NROW, NCOL), CHIP_DESIGN,
                                            dyes=0.2)

    # Sorting
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, dyes=(500, 200))
    assert_equal((200, 500), cgh_src.dyes)

    ## Bad tuple
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, dyes=(-500, 200))

    # Truncated
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, dyes=(50.30, 2.5))
    assert_equal((2, 50), cgh_src.dyes)

def test_parameters_outliers():
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, outliers=(0.3, 1.0))
    assert_equal((0.3, 1.0), cgh_src.outliers)

    # Default
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN)
    assert_equal((0.2, 0.5), cgh_src.outliers)

    # Not tuple
    assert_raises(TypeError, ArrayCGHSynth, (NROW, NCOL), CHIP_DESIGN,
                                            outliers=0.2)

    # Sorting
    cgh_src = ArrayCGHSynth((NROW, NCOL), CHIP_DESIGN, outliers=(0.5, 0.2))
    assert_equal((0.2, 0.5), cgh_src.outliers)

    # Bad tuple
    assert_raises(ValueError, ArrayCGHSynth,
                  (NROW, NCOL), CHIP_DESIGN, outliers=(-0.5, 0.2))

def test_created_signals():
    cgh_src = ArrayCGHSynth(geometry=(NROW, NCOL), design=CHIP_DESIGN)
    acgh = cgh_src.draw()

    for v in ('mask', 'test_signal', 'reference_signal',
              'start_base', 'end_base', 'chromosome', 'id',
              'true_test_signal'):
        assert_equal(acgh[v][~acgh['mask']], acgh.F[v])

    assert_equal(2., acgh.F['true_test_signal'])

def test_alterated_signal():
    cgh_src = ArrayCGHSynth(geometry=(NROW, NCOL), design=CHIP_DESIGN,
                            alterations={'1p12':[(4, 1.0)]},
                            cytostructure=CytoStructure(StringIO(CytoFileContent)))
    acgh = cgh_src.draw()
    acgh.sort()

    # Alterated
    assert_equal('P01', acgh.F['id'][0])
    assert_equal(4., acgh.F['true_test_signal'][0])

    # Normal
    assert_equal(2., acgh.F['true_test_signal'][1:])

