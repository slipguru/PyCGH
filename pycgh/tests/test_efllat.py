import numpy as np
from numpy.testing import *

# CGHDL utils -----------------------------------------------------------------
from ..analysis.efllat import (prox_squared_l1_bycol,
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
from ..analysis.efllat import prox_psi, prox_phi

def setup():
    return {'B': np.ones((100, 3)),
            'Theta': np.ones((3, 10)),
            'Y': np.ones((100, 10)),

            'eta': 1e-1,
            'zeta': 1e-1,
            'muw': np.ones((99, 1)),
            'lambda_': 1e-1,
            'tau': 1e-1,
            'bound': 1.0,
            'eps': 1e-3,
            'maxN': 1e3,
            'init': None}

def test_psi():
    params = prox_psi.func_code.co_varnames[:prox_psi.func_code.co_argcount]
    Zeta, gap, dvars = prox_psi(*(setup()[p] for p in params))

    assert_almost_equal(100.334, Zeta.sum(), 3)
    assert_almost_equal(0.034, gap, 3)

def test_phi():
    params = prox_phi.func_code.co_varnames[:prox_phi.func_code.co_argcount]
    Gamma, gap, dvars = prox_phi(*(setup()[p] for p in params))

    assert_almost_equal(10.631, Gamma.sum(), 3)
    assert_almost_equal(0.126, gap, 3)

# CGHDL main  algorithm -------------------------------------------------------
from ..analysis import cghDL

def setup_efllat():
    return {'Y': np.ones((100, 10)),

            'mu': 1.0,
            'tvw': np.ones((99, 1)),
            'lambda_': 1e-1,
            'tau': 1e-1,
            'J': 3,
            'theta_bound': 1.0,
            'eps': 1e-3,
            'maxN': 100,
            'maxK': 200,
            'initB': 'pca',
            'initTheta': None}

def test_efllat():
    params = cghDL.func_code.co_varnames[:cghDL.func_code.co_argcount]

    out = cghDL(*(setup_efllat()[p] for p in params))

    assert_almost_equal(10.0, out['Theta'].sum(), 3)
    assert_almost_equal(13.212, out['B'].sum(), 3)
    assert_almost_equal(4, out['conv'], 3)

    assert_almost_equal(0.124, out['gap_psi'], 3)
    assert_almost_equal(0.074, out['gap_phi'], 3)

# CGHDL parameter selection ---------------------------------------------------
from ..analysis import cghDL_BIC

from ..analysis.efllat import atoms_jumps
def test_atoms_jumps():
                #  1       2      2
    B = np.array([[1,      2,     3],
                  [1,      2.11,  1],
                  [1.01,   2,     3]])

    assert_equal(5, atoms_jumps(B, eps=1e-1))

def setup_efllat_bic():
    # Order J, lambda_, mu, tau (sorted)
    import itertools as it
    expected_params = it.product([1, 3], [1e-2, 1e-1], [1.0, 2.0], [1e-2, 1e-1])

    def cb(result, J, lambda_, mu, tau, BIC):
        # Check correct traversing rder
        assert_almost_equal(expected_params.next(), (J, lambda_, mu, tau))

    return {'Y': np.ones((100, 10)),

            'mu_range': [1.0, 2.0],
            'lambda_range': [1e-1, 1e-2],
            'tau_range': [1e-1, 1e-2],
            'J_range': [3, 1],
            'tvw': np.ones((99, 1)),
            'theta_bound': 1.0,
            'eps': 1e-3,
            'maxN': 100,
            'maxK': 200,
            'eps_jumps':1e-3,
            'initB': 'pca',
            'callback': cb}

def test_efllat_bic():
    params = cghDL_BIC.func_code.co_varnames[:cghDL_BIC.func_code.co_argcount]
    out = cghDL_BIC(*(setup_cghdl_bic()[p] for p in params))

    assert_almost_equal(-1817.533, out['BIC'], 3)
    assert_almost_equal(1e-2, out['lambda'], 3)
    assert_almost_equal(2.0, out['mu'], 3)
    assert_almost_equal(1e-2, out['tau'], 3)
