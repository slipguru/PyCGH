# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import itertools as it

import numpy as np

## Prox and Proj operators ----------------------------------------------------
from prox import prox_squared_l1_bycol # C-implementation
from utils import positive_box_projection, discrete_derivate, \
                  discrete_derivate_conj, interval_projection


## Auxiliary functions --------------------------------------------------------
def c_phi(Y, Theta0, B0, tau, eta, theta_bound):
    PTheta0 = np.empty_like(Theta0)
    positive_box_projection(Theta0, theta_bound, PTheta0)
    gap0_phi = (
        (0.5 * np.sum((Y - np.dot(B0, PTheta0))**2)) +
        (tau * (np.sum(np.sum(np.abs(PTheta0), axis=0)**2))) +
        (np.sum(PTheta0*PTheta0) / (2*eta)) +
        (np.sum(Theta0*Theta0) / (2*eta))
    )
    C_phi = np.sqrt(gap0_phi * 2.*eta)
    return C_phi

def c_psi(Y, Theta0, B0, lambda_, muw, zeta):
    L, J = B0.shape
    DB0 = np.empty((L-1, J))
    discrete_derivate(B0, DB0)
    gap0_psi = (
        (0.5 * np.sum((Y - np.dot(B0, Theta0))**2)) +
        (lambda_ * (np.sum(np.sum(np.abs(B0), axis=0)**2))) +
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
    #dual_prev = dual = np.inf
    gap = None

    for n in xrange(int(maxN)):
        #dual_prev = dual
        t_prev = t
        V1_prev, V2_prev, V3_prev = V1.copy(), V2.copy(), V3.copy()
        discrete_derivate_conj(U3, Dconj_U3)
        Zeta_aux = B - zeta*(np.dot(U1, Theta.T) + U2 + Dconj_U3)

        # Data fit:
        V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(Zeta_aux, Theta) - Y))

        # L1^2 norm
        grad = np.asfortranarray(U2 + gamma*Zeta_aux)
        prox_squared_l1_bycol(grad/gamma, V2, lambda_/gamma)
        V2 = grad - gamma*V2

        # Weighted Total variation
        discrete_derivate(Zeta_aux, DZeta)
        interval_projection(U3 + gamma*DZeta, mw, V3)

        # Solution Update
        discrete_derivate_conj(V3, Dconj_V3)
        Zeta = B - zeta*(np.dot(V1, Theta.T) + V2 + Dconj_V3)

        if not (n % 10):
            discrete_derivate(Zeta, DZeta)
            primal = (
                (0.5 * np.sum((Y - np.dot(Zeta, Theta))**2)) +            # Data fit
                (lambda_ * (np.sum(np.sum(np.abs(Zeta), axis=0)**2)) ) +  # L1^2
                (np.sum( muw * np.abs(DZeta)) ) +                         # TV
                (np.sum((Zeta - B)**2) / (2.*zeta))                       # Prox
            )
            dual = (
              ((np.sum(Zeta*Zeta) - BNorm2) / (2.*zeta)) +       # Prox*
              (0.5*np.sum(V1**2) + np.sum(V1*Y)) +               # Fit*
              ( np.sum(np.max(V2**2, axis=0))/(4.*lambda_))      # L1^2*
              # Tv not included because its dual is 0
            )

            gap = primal+dual

        # Restart!
        #if (dual > dual_prev):
        #    t = 1.
        #    U1 = V1.copy()
        #    U2 = V2.copy()
        #    U3 = V3.copy()
        #else:
        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        if gap <= (eps*eps)/(2.*zeta):
            break

    return Zeta, gap, (V1, V2, V3), n


def prox_phi(Theta, eta, B, Y, tau, bound, eps, maxN=1e5, init=None):
    """ Fixed B, Theta posbox"""

    J, S = Theta.shape
    L = B.shape[0]

    UBOUND = bound

    # Initializations
    ThetaNorm2 = np.sum(Theta*Theta)
    Gamma = np.empty_like(Theta)
    Gamma_aux = np.empty_like(Theta)
    PGamma = np.empty_like(Theta)

    if init is None:
        V1, V2, V3 = (np.zeros((L, S)), np.zeros((J, S)), np.zeros((J, S)))
    else:
        V1, V2, V3 = (V.copy() for V in init)
    U1, U2, U3 = V1.copy(), V2.copy(), V3.copy()
    V1_prev, V2_prev, V3_prev = (np.empty_like(V1), np.empty_like(V2),
                                 np.empty_like(V3))

    gamma = 1.0/(eta*(np.linalg.norm(np.dot(B.T, B)) + 2.0))
    t = 1.

    # Duality gap value
    #dual_prev = dual = np.inf
    gap = None

    for n in xrange(int(maxN)):
        #dual_prev = dual
        t_prev = t
        V1_prev, V2_prev, V3_prev = V1.copy(), V2.copy(), V3.copy()
        Gamma_aux = Theta - eta*(np.dot(B.T, U1) + U2 + U3)

        # Data fit
        V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(B, Gamma_aux) - Y))

        # Hard constraint (positivity and bound)
        grad = U2 + gamma*Gamma_aux
        positive_box_projection(grad, gamma*UBOUND, V2)
        V2 = grad - V2

        # L1^2 norm
        grad = np.asfortranarray(U3 + gamma*Gamma_aux)
        prox_squared_l1_bycol(grad/gamma, V3, tau/gamma)
        V3 = grad - gamma*V3

        # Solution Update
        Gamma = Theta - eta*(np.dot(B.T, V1) + V2 + V3)
        positive_box_projection(Gamma, UBOUND, PGamma)      # Pos+Box

        if not (n % 10):
            primal = (
                (0.5 * np.sum((Y - np.dot(B, PGamma))**2)) +           # Data Fit
                (tau * (np.sum(np.sum(np.abs(PGamma), axis=0)**2))) +  # L1^2
                (np.sum((PGamma - Theta)**2) / (2.*eta))               # Prox
                # Projection not included because PGamma is feasible
            )
            dual = (
              ((np.sum(Gamma*Gamma) - ThetaNorm2) / (2.*eta)) +       # Prox*
              (0.5*np.sum(V1**2) + np.sum(V1*Y)) +                    # Fit*
              (np.sum(np.max(V3**2, axis=0))/(4.*tau)) +              # L1^2*
              UBOUND * np.sum(np.clip(V2, 0., np.inf))                # ProjPos+Box*
            )

            gap = primal+dual

        ## Restart!
        #if (dual > dual_prev):
        #    t = 1.
        #    U1 = V1.copy()
        #    U2 = V2.copy()
        #    U3 = V3.copy()
        #else:
        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)
    
        if gap <= (eps*eps)/(2.*eta):
            break

    return PGamma, gap, (V1, V2, V3), n
