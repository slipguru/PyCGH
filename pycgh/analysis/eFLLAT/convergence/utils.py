# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Basic tools (TV and projections) """

import numpy as np

## Dictionary initialization options ------------------------------------------
def initB(Y, J, init_method='pca'):
    L, S = Y.shape

    if init_method in 'pca':
        TMP = 1./np.sqrt(S-1) * Y.T
        TMP -= np.mean(TMP, axis=0)
        U, d, Vt = np.linalg.svd(TMP, full_matrices=False)
        V = np.array(Vt.T)
        B0 = np.array(V[:, :J])
    elif init_method in 'rand':
        sampling = np.arange(S)
        np.random.shuffle(sampling)
        B0 = np.array(Y[:,sampling[:J]])

    return B0

## TV tools -------------------------------------------------------------------
def discrete_derivate(X, Y):
    X = np.asarray(X)
    Y = np.asarray(Y)
    Y[:] = X[1:] - X[:-1]
    return Y

def discrete_derivate_conj(X, Y):
    X = np.asarray(X)
    Y = np.asarray(Y)

    Y[0] = -X[0]
    Y[-1] = X[-1]
    Y[1:-1,:] = X[:-1] - X[1:]
    return Y

## Simple projections ---------------------------------------------------------
def interval_projection(X, w, Y):
    """
    Project x_i into the interval [-w_i, w_i] for each i.
    """
    X = np.asarray(X)
    np.clip(X, -w, w, out=Y)

def positive_box_projection(X, bound, Y):
    X = np.asarray(X)
    np.clip(X, 0.0, bound, out=Y)
    
def ball_proj(x, r=1.0):
    y = np.array(x, dtype=float) # copy

    y_norm = np.linalg.norm(y)
    if y_norm > r:
        y *= (r / y_norm)

    return y

def ball_proj_by_row(X, r=1.0):
    return np.apply_along_axis(ball_proj, 1, X, r)
    
## Prox l1 functions ----------------------------------------------------------
def prox_l1(x, t):
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.clip(np.abs(x) - t, a_min=0.0, a_max=np.inf)
