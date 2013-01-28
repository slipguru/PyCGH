import numpy as np
from numpy.testing import *

# CGHDL utils -----------------------------------------------------------------
from ..analysis.cghdl import pos_proj, ball_proj, symplex_proj

def test_pos_proj():
    assert_equal([1.0, 1.1, 0.0], pos_proj([1.0, 1.1, -0.1]))
    assert_equal([0.0, 0.0, 0.0], pos_proj([-1.0, -1.1, -0.1]))
    assert_equal([1.0, 1.1, 0.1], pos_proj([1.0, 1.1, 0.1]))

    assert_equal([1.0, 1.0, 0.0], pos_proj([1.0, 1.1, -0.1], r=1))
    assert_equal([0.0, 0.0, 0.0], pos_proj([-1.0, -1.1, -0.1], r=1))
    assert_equal([1.0, 1.0, 0.1], pos_proj([1.0, 1.1, 0.1], r=1))

    assert_equal([1.0, 1.1, 0.1], pos_proj([1.0, 1.1, 0.1], r=2))

def test_ball_proj():
    assert_equal([1.0, 0.0, 0.0], ball_proj([1.0, 0.0, 0.0]))
    assert_equal([1.0/np.sqrt(2.), 1.0/np.sqrt(2.), 0.0],
                  ball_proj([1.0, 1.0, 0.0]))
    assert_equal([1.0/np.sqrt(3.), 1.0/np.sqrt(3.), -1.0/np.sqrt(3)],
                  ball_proj([1.0, 1.0, -1.0]))

    assert_equal([1.0/(2*np.sqrt(1)), 0.0, 0.0],
                 ball_proj([1.0, 0.0, 0.0], r=1/2.))

def test_symplex_proj():
    assert_almost_equal(1, symplex_proj([0.3, 0.3, 0.3]).sum())
    assert_almost_equal(1, symplex_proj([0.3, 0.3, 0.5]).sum())
    assert_equal([0.5, 0.5, 0.0], symplex_proj([0.3, 0.3, -0.5]))


# CGHDL min  algorithm --------------------------------------------------------
from ..analysis import CGHDL

def test_default():
    cghDL = CGHDL()
