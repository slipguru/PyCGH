import numpy as np

## Prox and Proj operators ----------------------------------------------------
from utils import discrete_derivate, discrete_derivate_conj, \
                  interval_projection, prox_l1, ball_proj_by_row

## Auxiliary functions --------------------------------------------------------
def c_psi(Y, Theta0, B0, lambda_, mu, zeta):
    L, J = B0.shape
    DB0 = np.empty((L-1, J))
    discrete_derivate(B0, DB0)
    gap0_psi = (
        (0.5 * np.sum((Y - np.dot(B0, Theta0))**2)) +
        (lambda_ * np.sum(np.abs(B0))) +
        (np.sum( mu * np.abs(DB0))) +
        (np.sum(B0*B0) / zeta)
    )
    C_psi = np.sqrt(gap0_psi * 2.*zeta)
    return C_psi


def c_phi(Y, Theta0, B0, eta): # TAU ignored
    PTheta0 = np.empty_like(Theta0)
    PTheta0 = ball_proj_by_row(Theta0, 1.0)
    gap0_phi = (
        (0.5 * np.sum((Y - np.dot(B0, PTheta0))**2)) +
        (np.sum(PTheta0*PTheta0) / (2*eta)) +
        (np.sum(Theta0*Theta0) / (2*eta))
    )
    C_phi = np.sqrt(gap0_phi * 2.*eta)
    return C_phi


## FLLat prox steps -----------------------------------------------------------
def prox_psi(B, zeta, Theta, Y, mu, lambda_, eps, maxN=1e5, init=None):
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

    mw = (np.ones(L-1)*mu).reshape(-1, 1)

    if init is None:
        V1, V2, V3 = (np.zeros((L, S)), np.zeros((L, J)), np.zeros((L-1, J)))
    else:
        V1, V2, V3 = (V.copy() for V in init)
    U1, U2, U3 = V1.copy(), V2.copy(), V3.copy()
    V1_prev, V2_prev, V3_prev = (np.empty_like(V1), np.empty_like(V2),
                                 np.empty_like(V3))

    gamma = 1.0/(zeta * (np.linalg.norm(np.dot(Theta, Theta.T)) + 5.0))
    t = 1.
    # ISTA
    #gamma = 2.0/(zeta * (np.linalg.norm(np.dot(Theta, Theta.T)) + 5.0))

    # Duality gap value
    gaps_v1 = []
    gaps_v2 = []

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

        discrete_derivate(Zeta, DZeta)
        
        ### VERSION 1 #########################################################
        primal = (
            (0.5 * np.sum((Y - np.dot(Zeta, Theta))**2)) +     # Data fit
            (lambda_ * np.sum(np.abs(Zeta)) ) +                # L1               
            (np.sum( mu * np.abs(DZeta)) ) +                   # TV
            (np.sum((Zeta - B)**2) / (2.*zeta))                # Prox
        )
        dual = (
          ((np.sum(Zeta*Zeta) - BNorm2) / (2.*zeta)) +         # Prox*
          (0.5*np.sum(V1**2) + np.sum(V1*Y))                   # Fit*
          # Tv and L1 not included because their duals are 0
        )
        gaps_v1.append(primal+dual)
        
        ### VERSION 2 #########################################################
        gap = (
            (0.5 * np.sum((Y - np.dot(Zeta, Theta))**2)) +
            (lambda_ * np.sum(np.abs(Zeta)) ) +           
            (np.sum( mu * np.abs(DZeta)) )    +            
            (0.5*np.sum(V1**2) + np.sum(V1*Y)) +
            (np.sum((Zeta - B) * Zeta) / zeta) 
            # Tv and L1 not included because their duals are 0
        )        
        gaps_v2.append(gap)
        
        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)
        
        # ISTA
        #U1 = V1.copy()
        #U2 = V2.copy()
        #U3 = V3.copy()
        
        #if gap <= (eps*eps)/(2.*zeta):
            #break

    return Zeta, (gaps_v1, gaps_v2), (V1, V2, V3)
    
    
def prox_phi(Theta, eta, B, Y, eps, maxN=1e5, init=None):
    """ Fixed B, Theta l2"""

    J, S = Theta.shape
    L = B.shape[0]

    UBOUND = 1.0

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
    # ISTA
    #gamma = 2.0/(eta*(np.linalg.norm(np.dot(B.T, B)) + 1.0))

    # Duality gap value
    gaps_v1 = []
    gaps_v2 = []

    for n in xrange(int(maxN)):
        t_prev = t
        V1_prev, V2_prev = V1.copy(), V2.copy()
        Gamma_aux = Theta - eta*(np.dot(B.T, U1) + U2)

        # Data fit
        V1 = (1./(1. + gamma)) * (U1 + gamma*(np.dot(B, Gamma_aux) - Y))

        # Hard constraint (l2 bound)
        grad = U2 + gamma*Gamma_aux
        V2 = grad - ball_proj_by_row(grad, gamma*UBOUND)
        
        # Solution Update
        Gamma = Theta - eta*(np.dot(B.T, V1) + V2)
        PGamma = ball_proj_by_row(Gamma, UBOUND)      # l2 bound
        
        ### VERSION 1 #########################################################
        primal = (
            (0.5 * np.sum((Y - np.dot(B, PGamma))**2)) +           # Data Fit
            (np.sum((PGamma - Theta)**2) / (2.*eta))               # Prox
            # Projection not included because PGamma is feasible
        )
        dual = (
          ((np.sum(Gamma*Gamma) - ThetaNorm2) / (2.*eta)) +       # Prox*
          (0.5*np.sum(V1**2) + np.sum(V1*Y)) +                    # Fit*
          UBOUND * np.sum(np.sqrt(np.sum(V2**2, axis=1)))         # ProjL2*
        )
        gaps_v1.append(primal+dual)
        
        ### VERSION 2 #########################################################
        gap = (
            (0.5 * np.sum((Y - np.dot(B, PGamma))**2)) +
            (0.5*np.sum(V1**2) + np.sum(V1*Y)) +
            UBOUND * np.sum(np.sqrt(np.sum(V2**2, axis=1))) +
            (np.sum((PGamma - Theta) * PGamma) / eta) +
            ((np.sum(Gamma*Gamma) - np.sum(PGamma*PGamma)) / (2.*eta))
        )
        gaps_v2.append(gap)

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        
        # ISTA
        #U1 = V1.copy()
        #U2 = V2.copy()
        
        #if gap <= (eps*eps)/(2.*eta):
            #break
        
    return PGamma, (gaps_v1, gaps_v2), (V1, V2)