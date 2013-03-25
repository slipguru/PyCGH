import numpy as np
from numpy.testing import *

from ..datatypes import ArrayCGH
from ..analysis import cghnormaliter
from ..analysis.cghnormaliter import average_duplication

from test_arrayCGH import input

def test_input_data():
    aCGH = ArrayCGH(*input)
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)
    assert_equal(9, len(aCGH.names))

def test_average_no_duplication():
    aCGH = ArrayCGH(*input)
    aCGH = average_duplication(aCGH)
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)

def test_average_no_duplication_mask():
    mask = np.zeros(len(input[0]), dtype=bool)
    mask[:10] = True

    aCGH = ArrayCGH(*input, mask=mask)
    assert_equal(100, len(aCGH))
    assert_equal(90, aCGH.size)

    aCGH = average_duplication(aCGH)
    assert_equal(90, len(aCGH))
    assert_equal(90, aCGH.size)

def test_average_duplication():
    # Injecting duplication
    ids = input[0][:] ## Copy
    ids[:10] = ids[10:20]

    aCGH = ArrayCGH(ids, *input[1:])
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)

    aCGH = average_duplication(aCGH)
    assert_equal(90, len(aCGH))
    assert_equal(90, aCGH.size)

def test_average_duplication_mask():
    mask = np.zeros(len(input[0]), dtype=bool)
    mask[5:10] = True

    # Injecting duplication
    ids = input[0][:] ## Copy
    ids[:10] = ids[10:20]

    aCGH = ArrayCGH(ids, *input[1:], mask=mask)
    assert_equal(100, len(aCGH))
    assert_equal(95, aCGH.size)  # 5 are overlapping with masked

    aCGH = average_duplication(aCGH)
    assert_equal(90, len(aCGH))
    assert_equal(90, aCGH.size)

@dec.slow
def test_complete_run():
    aCGH = ArrayCGH(*input)         # No dup, no mask
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)

    aCGH = cghnormaliter(aCGH)
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)

    assert_equal(100, len(aCGH['cghnormaliter_ratio']))
    assert_equal(100, len(aCGH['cghnormaliter_call']))

@dec.slow
def test_complete_run_dup():
    # Injecting duplication
    ids = input[0][:] ## Copy
    ids[:10] = ids[10:20]

    aCGH = ArrayCGH(ids, *input[1:])  # Dup
    assert_equal(100, len(aCGH))
    assert_equal(100, aCGH.size)

    aCGH = cghnormaliter(aCGH)
    assert_equal(90, len(aCGH))
    assert_equal(90, aCGH.size)

    assert_equal(90, len(aCGH['cghnormaliter_ratio']))
    assert_equal(90, len(aCGH['cghnormaliter_call']))
