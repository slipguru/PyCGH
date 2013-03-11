import numpy as np
from numpy.testing import *

# CGHDL utils -----------------------------------------------------------------
from ..analysis.cghdl import (prox_squared_l1_bycol,
                              discrete_derivate, discrete_derivate_conj)

def test_prox_squared_l1():
    X = np.array([[1.0], [0.0], [0.0]])
    Y = np.empty_like(X)

    prox_squared_l1_bycol(X, Y, t=0.5)
    assert_equal([[0.5], [0.0], [0.0]], Y)

    prox_squared_l1_bycol(X, Y, t=0.3)
    assert_equal([[0.625], [0.0], [0.0]], Y)

    X[2,0] =  10.0
    prox_squared_l1_bycol(X, Y, t=0.3)
    assert_equal([[0.0], [0.0], [6.25]], Y)

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

def setup():
    return {'B': np.ones((100, 3)),
            'Theta': np.ones((3, 10)),
            'Y': np.ones((100, 10)),

            'eta': 1e-1,
            'zeta': 1e-1,
            'muw': np.ones((99, 1)),
            'mu': 1.0,
            'tvw': np.ones((99, 1)),
            'lambda_': 1e-1,
            'tau': 1e-1,
            'J': 1,
            'bound': 1.0,
            'eps': 1e-3,
            'maxN': 1e3,
            'maxK': 11,
            'init': None,
            'initB': 'pca'}

def test_psi():
    params = prox_psi.func_code.co_varnames[:prox_psi.func_code.co_argcount]
    Zeta, gaps, primals, duals, dvars = prox_psi(*(setup()[p] for p in params))

    assert_almost_equal(100.334, Zeta.sum(), 3)
    assert_almost_equal(1071.030, sum(gaps), 3)

def test_phi():
    params = prox_phi.func_code.co_varnames[:prox_phi.func_code.co_argcount]
    Gamma, gaps, primals, duals, dvars = prox_phi(*(setup()[p] for p in params))

    assert_almost_equal(10.631, Gamma.sum(), 3)
    assert_almost_equal(13.303, sum(gaps), 3)

# CGHDL main  algorithm -------------------------------------------------------
from ..analysis import cghDL, CGHDL
from ..analysis.orig_cghdl import cghDL

def test_cghdl():
    params = cghDL.func_code.co_varnames[:cghDL.func_code.co_argcount]

    #print params

    out = cghDL(*(setup()[p] for p in params))

    assert_almost_equal(0.073, out['Theta'].sum(), 3)
    assert_almost_equal(0.128, out['B'].sum(), 3)
    assert_almost_equal(10, out['conv'], 3)

def test_default():
    cghDL = CGHDL()
