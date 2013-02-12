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
    X = np.array([[1, 2, 3]]).T
    OUT = np.array([[1, 1]]).T
    Y = np.empty_like(OUT)
    assert_equal(OUT, discrete_derivate(X, Y))

def test_discrete_derivate_conj():
    X = np.array([[1, 2, 3]]).T
    OUT = np.array([[-1, -1, -1, 3]]).T
    Y = np.empty_like(OUT)
    assert_equal(OUT, discrete_derivate_conj(X, Y))

# CGHDL prox functions --------------------------------------------------------
from ..analysis.cghdl import prox_psi, prox_phi

def test_psi():
    B = np.ones((100, 3))
    Theta = np.ones((3, 10))
    Y = np.ones((100, 10))

    zeta = 1e-1
    muw = np.ones((99, 1))
    lambda_ = 1e-1
    eps = 1e-3
    maxN = 1e3

    Zeta, gaps, primals, duals, (V1, V2, V3) = prox_psi(B, zeta, Theta, Y,
                                                        muw, lambda_, eps,
                                                        maxN=maxN, init=None)
    assert_almost_equal(127.823, Zeta.sum(), 3)
    assert_almost_equal(107456.144, sum(gaps), 3)

def test_phi():
    B = np.ones((100, 3))
    Theta = np.ones((3, 10))
    Y = np.ones((100, 10))

    eta = 1e-1
    tau = 1e-1
    bound = 1.0
    eps = 1e-3
    maxN = 1e3

    Gamma, gaps, primals, duals, (V1, V2, V3) = prox_phi(Theta, eta, B, Y,
                                                        tau, bound, eps,
                                                        maxN=maxN, init=None)
    assert_almost_equal(10.606, Gamma.sum(), 3)
    assert_almost_equal(87.461, sum(gaps), 3)

# CGHDL main  algorithm -------------------------------------------------------
from ..analysis import CGHDL

def test_default():
    cghDL = CGHDL()
