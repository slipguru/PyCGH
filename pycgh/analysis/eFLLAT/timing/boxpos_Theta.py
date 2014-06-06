# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

## Prox and Proj operators ----------------------------------------------------
from utils import positive_box_projection

# Unchanged functions ---------------------------------------------------------
from cghdl import c_psi, prox_psi

## Auxiliary functions --------------------------------------------------------
def c_phi(Y, Theta0, B0, tau_ignored, eta, theta_bound): # TAU ignored
    PTheta0 = np.empty_like(Theta0)
    positive_box_projection(Theta0, theta_bound, PTheta0)
    gap0_phi = (
        (0.5 * np.sum((Y - np.dot(B0, PTheta0))**2)) +
        (np.sum(PTheta0*PTheta0) / (2*eta)) +
        (np.sum(Theta0*Theta0) / (2*eta))
    )
    C_phi = np.sqrt(gap0_phi * 2.*eta)
    return C_phi

## CGH DL steps ---------------------------------------------------------------
def prox_phi(Theta, eta, B, Y, tau_ignored, bound, eps, maxN=1e5, init=None): # TAU ignored
    """ Fixed B, Theta only posbox"""

    J, S = Theta.shape
    L = B.shape[0]

    UBOUND = bound

    # Initializations
    ThetaNorm2 = np.sum(Theta*Theta)
    Gamma = np.empty_like(Theta)
    Gamma_aux = np.empty_like(Theta)
    PGamma = np.empty_like(Theta)

    if init is None:
        V1, V2 = (np.zeros((L, S)), np.zeros((J, S)))
    else:
        V1, V2 = (V.copy() for V in init)
    U1, U2 = V1.copy(), V2.copy()
    V1_prev, V2_prev = (np.empty_like(V1), np.empty_like(V2))

    gamma = 1.0/(eta*(np.linalg.norm(np.dot(B.T, B)) + 1.0))
    t = 1.

    # Duality gap value
    gap = None

    for n in xrange(int(maxN)):
        t_prev = t
        V1_prev, V2_prev = V1.copy(), V2.copy()
        Gamma_aux = Theta - eta*(np.dot(B.T, U1) + U2)

        # Data fit
        V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(B, Gamma_aux) - Y))

        # Hard constraint (positivity and bound)
        grad = U2 + gamma*Gamma_aux
        positive_box_projection(grad, gamma*UBOUND, V2)
        V2 = grad - V2
        
        # Solution Update
        Gamma = Theta - eta*(np.dot(B.T, V1) + V2)
        positive_box_projection(Gamma, UBOUND, PGamma)      # Pos+Box

        if not (n % 10):
            primal = (
                (0.5 * np.sum((Y - np.dot(B, PGamma))**2)) +           # Data Fit
                (np.sum((PGamma - Theta)**2) / (2.*eta))               # Prox
                # Projection not included because PGamma is feasible
            )
            dual = (
              ((np.sum(Gamma*Gamma) - ThetaNorm2) / (2.*eta)) +       # Prox*
              (0.5*np.sum(V1**2) + np.sum(V1*Y)) +                    # Fit*
              UBOUND * np.sum(np.clip(V2, 0., np.inf))                # ProjPos+Box*
            )

            gap = primal+dual

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)

        if gap <= (eps*eps)/(2.*eta):
            break

    return PGamma, gap, (V1, V2), n
