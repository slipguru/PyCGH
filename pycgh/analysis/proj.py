# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np
from numpy.testing import *

## Projections ----------------------------------------------------------------
def pos_proj(x, r=np.inf):
    return np.clip(x, a_min=0.0, a_max=r)

def ball_proj(x, r=1.0):
    y = np.array(x, dtype=float) # copy

    y_norm = np.linalg.norm(y)
    if y_norm > r:
        y *= (r / y_norm)

    return y

def pos_ball_proj(x, r=1.0):
    x = np.asarray(x)
    if np.alltrue(x >= 0):
        return ball_proj(x, r)
    elif np.alltrue(x < 0):
        return np.zeros_like(x)

    return ball_proj(x.clip(min=0.0), r) # copy

def simplex_proj(x, r=1.0):
    x = np.asarray(x)
    if x.sum() == r and np.alltrue(x >= 0):
        return np.array(x, dtype=float) # copy

    x_s = np.sort(x)[::-1]
    xtmp = (x_s.cumsum() - r) / np.arange(1, len(x) + 1)
    k = np.where(xtmp < x_s)[0][-1]

    return (x - np.max(xtmp[k], 0)).clip(min=0)

def prox_l1(x, t):
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.clip(np.abs(x) - t, a_min=0.0, a_max=np.inf)

## Testing --------------------------------------------------------------------
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
