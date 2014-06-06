# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

## Prox and Proj operators ----------------------------------------------------
from utils import discrete_derivate, discrete_derivate_conj, \
                  interval_projection, prox_l1

# Unchanged functions ---------------------------------------------------------
from l12_B import c_phi, prox_phi

## Auxiliary functions --------------------------------------------------------
def c_psi(Y, Theta0, B0, lambda_, muw, zeta):
    L, J = B0.shape
    DB0 = np.empty((L-1, J))
    discrete_derivate(B0, DB0)
    gap0_psi = (
        (0.5 * np.sum((Y - np.dot(B0, Theta0))**2)) +
        (lambda_ * np.sum(np.abs(B0))) +
        (np.sum( (muw) * np.abs(DB0))) +
        (np.sum(B0*B0) / zeta)
    )
    C_psi = np.sqrt(gap0_psi * 2.*zeta)
    return C_psi


## CGH DL steps ---------------------------------------------------------------
def prox_psi(B, zeta, Theta, Y, muw, lambda_, eps, maxN=1e5, init=None):
    """ Fixed Theta """

    L, J = B.shape
    S = Theta.shape[1]

    # Initializations
    BNorm2 = np.sum(B*B)
    Zeta = np.empty_like(B)
    Zeta_aux = np.empty_like(B)
    Dconj_U3 = np.empty_like(B)
    Dconj_V3 = np.empty_like(B)
    DZeta = np.empty((L-1, J))

    mw = muw.ravel().reshape(-1, 1)

    if init is None:
        V1, V2, V3 = (np.zeros((L, S)), np.zeros((L, J)), np.zeros((L-1, J)))
    else:
        V1, V2, V3 = (V.copy() for V in init)
    U1, U2, U3 = V1.copy(), V2.copy(), V3.copy()
    V1_prev, V2_prev, V3_prev = (np.empty_like(V1), np.empty_like(V2),
                                 np.empty_like(V3))

    gamma = 1.0/(zeta * (np.linalg.norm(np.dot(Theta, Theta.T)) + 5.0))
    t = 1.

    # Duality gap value
    gap = None

    for n in xrange(int(maxN)):
        t_prev = t
        V1_prev, V2_prev, V3_prev = V1.copy(), V2.copy(), V3.copy()
        discrete_derivate_conj(U3, Dconj_U3)
        Zeta_aux = B - zeta*(np.dot(U1, Theta.T) + U2 + Dconj_U3)

        # Data fit:
        V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(Zeta_aux, Theta) - Y))

        # L1 norm
        grad = U2 + gamma*Zeta_aux
        interval_projection(grad, lambda_, V2)

        # Weighted Total variation
        discrete_derivate(Zeta_aux, DZeta)
        interval_projection(U3 + gamma*DZeta, mw, V3)

        # Solution Update
        discrete_derivate_conj(V3, Dconj_V3)
        Zeta = B - zeta*(np.dot(V1, Theta.T) + V2 + Dconj_V3)

        if not (n % 10):
            discrete_derivate(Zeta, DZeta)
            primal = (
                (0.5 * np.sum((Y - np.dot(Zeta, Theta))**2)) +     # Data fit
                (lambda_ * np.sum(np.abs(Zeta)) ) +                # L1               
                (np.sum( muw * np.abs(DZeta)) ) +                  # TV
                (np.sum((Zeta - B)**2) / (2.*zeta))                # Prox
            )
            dual = (
              ((np.sum(Zeta*Zeta) - BNorm2) / (2.*zeta)) +         # Prox*
              (0.5*np.sum(V1**2) + np.sum(V1*Y))                   # Fit*
              # Tv and L1 not included because their duals are 0
            )

            gap = primal+dual

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        if gap <= (eps*eps)/(2.*zeta):
            break

    return Zeta, gap, (V1, V2, V3), n