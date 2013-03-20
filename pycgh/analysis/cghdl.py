# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import itertools as it

import numpy as np

## Prox operators -------------------------------------------------------------
from prox import prox_squared_l1_bycol # C-implementation

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
        discrete_derivate(Zeta, DZeta)

        if not (n % 10):
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

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        if gap <= (eps*eps)/(2.*zeta):
            break

    return Zeta, gap, (V1, V2, V3)


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
    gap = None

    for n in xrange(int(maxN)):
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

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1 = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2 = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3 = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        if gap <= (eps*eps)/(2.*eta):
            break

    return PGamma, gap, (V1, V2, V3)

### CGHDL MAIN FUNCTION #######################################################
def _initB(Y, J, init_method='pca'):
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

def cghDL(Y, J, lambda_, mu, tau, theta_bound=1.0, tvw=None,
          initB='pca', initTheta=None, maxK=200, maxN=100, eps=1e-3):

    L, S = Y.shape

    #### B Initialization
    try:
        assert initB.shape == (L, J)
        B0 = initB.copy()
    except AttributeError:
        assert initB in ['pca', 'rand']
        B0 = _initB(Y, J, initB)

    #### Theta Initialization
    theta_bound = float(theta_bound)
    try:
        assert initTheta.shape == (J, S)
        Theta0 = initTheta.copy()
    except AttributeError:
        Theta0 = np.ones((J, S)) * theta_bound

    #### Total variation weigths
    if tvw is None:
        w = np.ones((L-1, 1))
    else:
        w = tvw.reshape(L-1, 1)

    #### Precision
    p = 2.
    eta = 1./(S*J)
    zeta = 1./(L*J)

    # Parameters normalization
    tau *= (L/float(J))      # s, j -- /J
    lambda_ *= (S/float(J))  # l, j -- /S
    mu *= (S/float(J))       # l, j -- /S

    #### Starting Duality Gap for PHI (with Vi=0, B fixed and Gamma=Theta0)
    PTheta0 = np.empty_like(Theta0)
    positive_box_projection(Theta0, theta_bound, PTheta0)
    gap0_phi = (
        (0.5 * np.sum((Y - np.dot(B0, PTheta0))**2)) +
        (tau * (np.sum(np.sum(np.abs(PTheta0), axis=0)**2))) +
        (np.sum(PTheta0*PTheta0) / (2*eta)) +
        (np.sum(Theta0*Theta0) / (2*eta))
    )
    C_phi = np.sqrt(gap0_phi * 2.*eta)

    #### Starting Duality Gap for PSI (with Vi=0, Theta fixed and Zeta=B0)
    DB0 = np.empty((L-1, J))
    discrete_derivate(B0, DB0)
    gap0_psi = (
        (0.5 * np.sum((Y - np.dot(B0, Theta0))**2)) +
        (lambda_ * (np.sum(np.sum(np.abs(B0), axis=0)**2))) +
        (np.sum( (mu*w) * np.abs(DB0))) +
        (np.sum(B0*B0) / zeta)
    )
    C_psi = np.sqrt(gap0_psi * 2.*zeta)

    dual_var_phi = None
    dual_var_psi = None
    B, Theta = B0, Theta0

    B_diff = None
    Theta_diff = None
    B_prev = np.empty_like(B)
    Theta_prev = np.empty_like(Theta)

    for k in xrange(int(maxK)):
        B_prev = B.copy()
        Theta_prev = Theta.copy()

        epsk = 1. / ((k+1)**p)
        eta = 1. / ((k+1)**p)
        zeta = 1. / ((k+1)**p)

        # Theta update
        Theta, gap_phi, dual_var_phi = prox_phi(Theta, eta, B, Y, tau,
                                                bound=theta_bound,
                                                eps=C_phi*epsk,
                                                maxN=maxN,
                                                init=dual_var_phi)

        # B update
        B, gap_psi, dual_var_psi = prox_psi(B, zeta, Theta, Y, mu*w, lambda_,
                                            eps=C_psi*epsk,
                                            maxN=maxN,
                                            init=dual_var_psi)

        B_diff = np.sum((B - B_prev)**2)/np.sum(B_prev**2)
        Theta_diff = np.sum((Theta - Theta_prev)**2)/np.sum(Theta_prev**2)

        # convergence
        if (B_diff <= eps and Theta_diff <= eps):
            break

    return {'B': B, 'Theta': Theta, 'conv': k,
            'gap_phi': gap_phi, 'gap_psi': gap_psi}


### CGHDL BIC-BASED PARAMETER SEARCHING #######################################
def atoms_jumps(B, eps=1e-3):
    """ Counts the number of non-zero levels (up to eps tolerance). """
    ajumps = 0.0
    bs = np.empty_like(B[:,0])

    for j in xrange(B.shape[1]):
        bs = B[:,j].copy()
        bs.sort()

        ajumps += (np.diff(bs) > eps).sum() + 1.0

    return ajumps


def cghDL_BIC(Y, J_range, lambda_range, mu_range, tau_range,
              theta_bound=1.0, tvw=None, initB='pca',
              maxK=200, maxN=100, eps=1e-3, eps_jumps=1e-3,
              callback=None):

    L, S = Y.shape
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)

    # Parameters dimensions
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)

    # CGHDL params
    params = {'theta_bound':theta_bound, 'tvw':tvw,
              'maxK':maxK, 'maxN':maxN, 'eps':eps}

    # Initialization
    BIC_min = np.inf
    B_res, Theta_res = None, None
    J_res, mu_res, lambda_res, tau_res = None, None, None, None

    for J in J_range:
        # For each new J value we start from following B0 and Theta0
        B0 = _initB(Y, J)
        Theta0 = np.ones((J, S)) * theta_bound


        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):

            result = cghDL(Y, J,
                           lambda_=lambda_range[j],
                           mu=mu_range[i],
                           tau=tau_range[k],
                           initB=B0, initTheta=Theta0, **params)
            B = result['B']
            Theta = result['Theta']

            # BIC calculation
            fit = np.sum((Y - np.dot(B, Theta))**2.) / (S*L)
            jumpsB = atoms_jumps(B, eps_jumps)

                        # fit                 # complexity
            BIC = (S*L * np.log(fit)) + (jumpsB * np.log(S*L))

            if BIC < BIC_min:
                BIC_min = BIC
                B_res, Theta_res = B, Theta
                J_res, lambda_res, mu_res, tau_res = (J, lambda_range[j],
                                                         mu_range[i],
                                                         tau_range[k])
            # Callback execution
            if not callback is None:
                callback(result, J,
                         lambda_range[j], mu_range[i], tau_range[k], BIC)

    # Return the best result
    return {'B':B_res, 'Theta':Theta_res, 'J':J_res,
            'lambda':lambda_res, 'mu':mu_res, 'tau':tau_res, 'BIC':BIC_min}
