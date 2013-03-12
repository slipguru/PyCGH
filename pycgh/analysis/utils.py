import numpy as np

import scipy.lib.blas as blas

def apply_by_row(func, x, *args):
    return np.apply_along_axis(func, 1, x, *args)

def apply_by_col(func, x, *args):
    return np.apply_along_axis(func, 0, x, *args)

def _norm2(q):
    q = np.asarray(q)
    nrm2, = blas.get_blas_funcs(['nrm2'], [q])
    return nrm2(q)

def box_proj(x, r=1.0):
    x = np.asarray(x, dtype=float)

    x_proj = x.copy()
    x_proj[x < 0] = 0.0
    x_proj[x > r] = r

    return x_proj

def ball_proj(x, r=1.0):
    x = np.asarray(x, dtype=float)

    x_norm = np.linalg.norm(x)
    if x_norm <= r:
        return x
    else:
        return (r*x) / x_norm

def pos_proj(x):
    x = np.asarray(x, dtype=float)

    neg_x = (x < 0)
    x_proj = x.copy()
    x_proj[neg_x] = 0.0

    return x_proj

def neg_proj(x):
    x = np.asarray(x, dtype=float)

    pos_x = (x > 0)
    x_proj = x.copy()
    x_proj[pos_x] = 0.0

    return x_proj

def symplex_proj(x):
    # Projection onto the standar symplex of R^n

    n = len(x)

    x = np.asarray(x, dtype=float)
    #sorted_idx = np.argsort(x)[::-1] # sorted reversed
    #x_s = x[sorted_idx]
    x_s = np.sort(x)[::-1]

    xtmp = (x_s.cumsum()-1)/np.arange(1,n+1)

    k = np.max(np.where(xtmp < x_s)[0])

    t = np.maximum(xtmp[k],0)

    px = np.maximum(x - t, 0)

    return px

def symplex_suppfunc(u):

    u = np.asarray(u, dtype=float)

    up=np.maximum(u,0)

    return up.sum()



def pos_l2_proj(x, r=1.0):
    """ Projection on the positive orthant of the l2 ball of radius r.

    Parameters
    ----------
    x : array_like, shape (N, )
        Vector to project.
    r : floar
        L2 ball radius.

    Returns
    -------
    xp : array_like, shape (N, )
         Projected vector.

    """
    x = np.asarray(x, dtype=float)

    pos_x = (x > 0)
    neg_x = (x <= 0)

    # All negative components
    if neg_x.sum() == len(x):
        return np.zeros_like(x)

    # Positive orthant
    if pos_x.sum() == len(x):
        x_norm = np.linalg.norm(x)
        if x_norm <= r:
            return x
        else:
            return (r*x) / x_norm

    # Otherwise
    x_proj = np.empty_like(x)
    x_proj[neg_x] = 0.0

    x_pos_norm = np.linalg.norm(x[pos_x])
    if x_pos_norm <= r:
        x_proj[pos_x] = x[pos_x]
    else:
        x_proj[pos_x] = (r*x[pos_x]) / x_pos_norm

    return x_proj

def interval_projection(x, w):
    """
    Project x_i into the interval [-w_i, w_i] for each i.
    """
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    assert np.all(w >= 0)

    return np.clip(x, -w, w)

def _soft_thresholding(x, t):
    return x - interval_projection(x, t)

def soft_thresholding(x, t):
    assert t >=0
    return np.sign(x) * np.clip(np.abs(x) - t, 0.0, np.inf)

def prox_l1_squared_norm(x, lambda_):
    """
    Proximal operator for the squared l1 norm.

    The computed prox has the form:
        prox_{lambda ||.||_1^2}(x)

    Parameters
    ----------
    x : array_like, shape (N,)

    Returns
    -------
    x_prox : array_like, shape(N,)

    """
    x = np.asarray(x, dtype=float)
    sorted_idx = np.argsort(np.abs(x))[::-1] # sorted reversed
    x_s = np.abs(x[sorted_idx])

    right_cond = 2.*lambda_ * x_s.cumsum()
    left_cond = 1. + 2.*lambda_*np.arange(1, len(x)+1)

    l = 0
    for l in xrange(len(x)-1):
        if ( (left_cond[l] * x_s[l+1] <= right_cond[l]) and
             (left_cond[l] * x_s[l] > right_cond[l]) ):
            break

    rho = right_cond[l] / left_cond[l] # last computed conditions
    #print rho
    return soft_thresholding(x, rho)

#def prox_l1_squared_norm_bycol(X, lambda_):
#    X = np.asfortranarray(X, dtype=float)
#    #sorted_idx = np.argsort(np.abs(X), axis=0)[::-1]
#    X_s = np.sort(np.abs(X), axis=0)[::-1]
#
#    L, J = X.shape
#
#    right_cond = X_s.cumsum(axis=0)     # (L-1, J)
#    left_cond = (1./(2.*lambda_) + np.arange(1, X.shape[0]+1))  # L
#    left_cond.shape = (-1, 1)
#
#    last_conds = np.logical_and(left_cond[:-1] * X_s[1:] <= right_cond[:-1],
#                                left_cond[:-1] * X_s[:-1] > right_cond[:-1])
#
#    last_conds = np.r_[last_conds, np.ones((1, X.shape[1]), dtype=bool)]
#
#    last_conds = np.where(last_conds)
#    rho = right_cond[(last_conds[0][:J], last_conds[1][:J])]/left_cond[last_conds[0][:J]].T
#    return np.sign(X) * np.clip(np.abs(X) - rho, 0.0, np.inf)

def discr_derivate(x):
    x = np.asarray(x)
    return x[1:] - x[:-1]

def discr_derivate_conj(x):
    x = np.asarray(x)
    return np.r_[-x[0], x[:-1] - x[1:], x[-1]]
