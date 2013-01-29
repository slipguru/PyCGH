import numpy as np
from numpy.testing import *

# CGHDL utils -----------------------------------------------------------------
from ..analysis.cghdl import pos_proj, ball_proj, pos_ball_proj, simplex_proj
from ..analysis.cghdl import prox_l1, prox_squared_l1
from ..analysis.cghdl import discrete_derivate, discrete_derivate_conj

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

def test_pos_ball_proj():
    assert_equal([1.0, 0.0, 0.0], pos_ball_proj([1.0, 0.0, 0.0]))
    assert_equal([1.0, 0.0, 0.0], pos_ball_proj([1.0, -1.0, 0.0]))

    assert_equal([1.0/np.sqrt(2.), 1.0/np.sqrt(2.), 0.0],
                  pos_ball_proj([1.0, 1.0, 0.0]))
    assert_equal([1.0/np.sqrt(2.), 1.0/np.sqrt(2.), 0.0],
                  pos_ball_proj([1.0, 1.0, -1.0]))

def test_simplex_proj():
    assert_almost_equal(1, simplex_proj([0.3, 0.3, 0.3]).sum())
    assert_almost_equal(1, simplex_proj([20, 40, 10.]).sum())
    assert_almost_equal(1, simplex_proj([0.3, 0.3, 0.5]).sum())

    assert_almost_equal(2., simplex_proj([0.3, 0.3, -0.5], r=2.).sum())

def test_prox_l1():
    assert_equal([0.5, 0.0, 0.0], prox_l1([1.0, 0.0, 0.0], t=0.5))
    assert_equal([0.7, 0.0, 0.0], prox_l1([1.0, 0.0, 0.0], t=0.3))
    assert_equal([0.7, 0.0, 9.7], prox_l1([1.0, 0.0, 10.0], t=0.3))

def test_prox_squared_l1():
    assert_equal([0.5, 0.0, 0.0], prox_squared_l1([1.0, 0.0, 0.0], t=0.5))
    assert_equal([0.625, 0.0, 0.0], prox_squared_l1([1.0, 0.0, 0.0], t=0.3))
    assert_equal([0.625, 0.0, 9.625], prox_squared_l1([1.0, 0.0, 10.0], t=0.3))
    
def test_discrete_derivate():
    assert_equal([1, 1], discrete_derivate([1, 2, 3]))
    assert_equal([], discrete_derivate([1]))

def test_discrete_derivate_conj():
    assert_equal([-1, -1, -1, 3], discrete_derivate_conj([1, 2, 3]))
    assert_equal([-1, 1], discrete_derivate_conj([1]))

# CGHDL min  algorithm --------------------------------------------------------
from ..analysis import CGHDL

def test_default():
    cghDL = CGHDL()
