# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

# CGH DL utils ----------------------------------------------------------------
## Projections -----
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

## Prox operators -----
def prox_l1(x, t):
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.clip(np.abs(x) - t, a_min=0.0, a_max=np.inf)

from prox import prox_squared_l1 as psl1
def prox_squared_l1(x, t):
    return psl1(np.asarray(x), float(t))

## Derivate -----
def discrete_derivate(x):
    x = np.asarray(x)
    return x[1:] - x[:-1]

def discrete_derivate_conj(x):
    x = np.asarray(x)
    return np.r_[-x[0], x[:-1] - x[1:], x[-1]]
    
## Supp functions ------
#def symplex_suppfunc(u):
#    u = np.asarray(u, dtype=float)
#    up = np.maximum(u, 0)
#    return up.sum()

class CGHDL(object):
    def __init__(self):
        pass
