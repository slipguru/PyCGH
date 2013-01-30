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

## Derivate by column -----
#def discrete_derivate(x):
#    x = np.asarray(x)
#    return x[1:] - x[:-1]

def discrete_derivate(X, Y):
    X = np.asarray(X)
    Y = np.asarray(Y)
    Y[:] = X[1:] - X[:-1]

#def discrete_derivate_conj(x):
    #x = np.asarray(x)
    #return np.r_[-x[0], x[:-1] - x[1:], x[-1]]

def discrete_derivate_conj(X, Y):
    X = np.asarray(X)
    Y = np.asarray(Y)

    Y[0] = -X[0]
    Y[-1] = X[-1]
    Y[1:-1,:] = X[:-1] - X[1:]

## Supp functions ------
#def symplex_suppfunc(u):
#    u = np.asarray(u, dtype=float)
#    up = np.maximum(u, 0)
#    return up.sum()

################## LIMBO ######################################################
def apply_by_row(func, x, *args):
    return np.apply_along_axis(func, 1, x, *args)

from prox import prox_squared_l1_bycol

def apply_by_col(func, x, *args):
    return np.apply_along_axis(func, 0, x, *args)

def interval_projection(X, w, Y):
    """
    Project x_i into the interval [-w_i, w_i] for each i.
    """
    X = np.asarray(X)
    np.clip(X, -w, w, out=Y)

# CGH DL steps ----------------------------------------------------------------
from scipy.linalg.blas import get_blas_funcs
dgemm = get_blas_funcs('gemm')

def prox_psi(B, zeta, Theta, Y, muw, lambda_, eps, maxN=1e5, init=None):
    """ Fixed Theta """

    L, J = B.shape
    S = Theta.shape[1]

    # Initializations
    BNorm2 = np.sum(B*B)
    Zeta = np.empty_like(B, order='F')
    Zeta_aux = np.empty_like(B, order='F')
    Dconj_U3 = np.empty_like(B, order='F')
    Dconj_V3 = np.empty_like(B, order='F')
    DZeta = np.empty((L-1, J), order='F')
    mw = muw.ravel().reshape(-1, 1) # not diag...

    if init is None:
        V1, V2, V3 = (np.zeros((L, S), order='F'),
                      np.zeros((L, J), order='F'),
                      np.zeros((L-1, J), order='F'))
    else:
        V1, V2, V3 = init
    U1, U2, U3 = V1.copy(), V2.copy(), V3.copy()
    V1_prev, V2_prev, V3_prev = (np.empty_like(V1, order='F'),
                                 np.empty_like(V2, order='F'),
                                 np.empty_like(V3, order='F') )

    gamma = 1.0/(zeta * (np.linalg.norm(np.dot(Theta, Theta.T)) + 5.0))
    t = 1.

    # GAPS values
    gaps = list()
    primals = list()
    duals = list()

    for n in xrange(int(maxN)):
        t_prev = t
        V1_prev[:], V2_prev[:], V3_prev[:] = V1, V2, V3
        discrete_derivate_conj(U3, Dconj_U3)
        Zeta_aux[:] = B - zeta*(np.dot(U1, Theta.T) + U2 + Dconj_U3)

        # Data fit:
        #    V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(Zeta_aux, Theta) - Y))
        V1[:] = (U1 - gamma*Y)
        dgemm((gamma/(1. + gamma)), Zeta_aux, Theta,  (1./(1. + gamma)), V1, overwrite_c=True)

        # L1^2 norm
        grad = U2 + gamma*Zeta_aux
        prox_squared_l1_bycol(grad/gamma, V2, lambda_/gamma)
        V2[:] = grad - gamma * V2

        # Weighted Total variation
        discrete_derivate(Zeta_aux, DZeta)
        interval_projection(U3 + gamma*DZeta, mw, V3)

        # Solution Update
        discrete_derivate_conj(V3, Dconj_V3)
        #Zeta[:] = B - zeta*(np.dot(V1, Theta.T) + V2 + Dconj_V3)
        Zeta[:] = V2 + Dconj_V3 - B/zeta
        dgemm(-zeta, V1, Theta, -zeta, Zeta, trans_b=True, overwrite_c=True)
        discrete_derivate(Zeta, DZeta)

        if n%10==0:
          primal = (
            (0.5 * np.sum((Y - np.dot(Zeta, Theta))**2)) +            # Data fit
            (lambda_ * (np.sum(np.sum(np.abs(Zeta), axis=0)**2)) ) +  # L1^2
            #(np.sum(np.dot(mw, np.abs(DZeta)))) +                     # TV
            (np.sum( muw * np.abs(DZeta)) ) +
            (np.sum((Zeta - B)**2) / (2.*zeta))                       # Prox
          )
          dual = (
            ((np.sum(Zeta*Zeta) - BNorm2) / (2.*zeta)) +       # Prox*
            (0.5*np.sum(V1**2) + np.sum(V1*Y)) +               # Fit*
            ( np.sum(np.max(V2**2, axis=0))/(4.*lambda_))      # L1^2*
            # Tv not included because its dual is 0
          )

          gap = primal+dual
          gaps.append(gap)
          primals.append(primal)
          duals.append(dual)

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1[:] = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2[:] = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3[:] = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        #if (n == 0) or ((n+1) % 1000 == 0):
        #    print 'Prox Psi #it: %d' % (n+1)
        #    print '   Gap: %.5f (th %.5e)' % (gap, (eps*eps)/(2.*zeta))
        #    print '   Primal: %.5f' % primal
        #    print '   Dual: %.5f' % dual

        if gap <= (eps*eps)/(2.*zeta):
            #print '##-Exit with gap %.5e (#%d it.)-##' % (gap, (n+1))
            break

    #if n == (maxN-1):
    #    print '##-Exit with gap %.5e (reached maximum #it %d)-##' % (gap, maxN)

    return Zeta, gaps, primals, duals, (V1, V2, V3)

class CGHDL(object):
    def __init__(self):
        pass
