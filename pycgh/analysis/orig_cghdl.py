##########################################################################
# First Draft implemetation of the DL-based method for aCGH segmentation #
##########################################################################

import numpy as np

from scipy import signal

from utils import *

import pylab as pl


# S aCGH data are assumed stacked on the columns of the LxS-dimensional matrix
#   Y = [y_1  ... y_S], for each s=1,...,S y_s in R^L
# and composed by the linear composition of J basic patterns (called atoms),
# stacked on the columns of the LxJ-dimensional matrix
#   B = [beta_1 ... beta_J], for each j=1,...,J beta_j in R^L.
#
# Our objective is to reconstruct each sample by a sparse combination
# of atoms
#   Y = B*Theta + E,
# where E is an LxS -dimensional noise matrix and the JxS-dimensional matrix
# Theta contains the weights associated with the reconstruction
#   Theta = [theta_1 ... theta_S], for each s=1,...,S theta_s in R^J
# such that
#   y_s = sum(beta_j * theta_s[j])_{j=1,...,J} + eps

def get_data(J=30, S=100, L=1000):
    indexes = np.array_split(np.arange(L), J)

    # Testing Data
    B = np.empty((L, J))
    for j, idxs in enumerate(indexes):
        profile = np.zeros(L)
        profile[idxs] = 1. * (-1 if j%2 else 1)

        B[:, j] = profile

    Theta = np.zeros((J, S))
    for j, idxs in enumerate(np.array_split(np.arange(S), J)):
        Theta[j, idxs] = 10.

    #for s in xrange(S):
        #idxs = np.random.randint(J, size=3) ### if replaced less atom used!
        #Theta[s%J, s] = 10.

    Y = np.dot(B, Theta) + np.random.random((L, S))

    return Y, Theta, B

#def proj_grad(Theta, eta, B, Y, maxN=1e6):
#    X = np.zeros_like(Theta)
#    Gamma = np.zeros_like(Theta)
#
#    L, J = B.shape
#
#    BTB = np.dot(B.T, B)
#    BTY = np.dot(B.T, Y)
#    I = np.eye(J)
#    gamma_n = 1.0 / (eta*np.linalg.norm(BTB)**2 + 1.0)
#
#    t = t_prev = 1.
#    X_prev = X.copy()
#    func_values = list()
#    n = 0
#    while n < maxN: ###<--
#        n += 1
#
#        grad = np.dot(eta*BTB + I, Gamma) - (eta*BTY + Theta)
#        arg = Gamma - gamma_n*grad
#        X = np.apply_along_axis(pos_proj, 1, arg)#, RADIUS) ##########
#        #X = arg # No Proj
#
#        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
#        Gamma = X + ((t_prev - 1.) / t) * (X - X_prev)
#
#        t_prev = t
#        X_prev = X.copy()
#
#        f_value = (( 0.5 * np.linalg.norm(Y - np.dot(B, X))**2 ) +
#                   ( (0.5/eta) * np.linalg.norm(X - Theta)**2))
#        func_values.append(f_value)
#
#        if n % 1000 == 1:
#            print n #, gap, (eps*eps)/(2.*eta)
#
#    return X, func_values

def prox_phi(Theta, eta, B, Y, tau, bound, eps, maxN=1e5, init=None,
             constr='posbox'):
             #eps, maxN=1e5, init=None, constr='posl2',
             #bound=1.0):
    """ Fixed B """

    J, S = Theta.shape
    L = B.shape[0]

    # Input check (not optimized!)
    assert constr in ['pos', 'posl2', 'posbox', 'l2', 'none']

    if constr == 'posl2':
        RADIUS = bound#1.#S*10.
    elif constr == 'posbox':
        UBOUND = bound#s10.#0.

    # Initializations
    ThetaNorm2 = np.sum(Theta*Theta)
    Gamma = np.empty_like(Theta, order='C')
    Gamma_aux = np.empty_like(Theta, order='C')
    PGamma = np.empty_like(Theta, order='C')   

    if init is None:
        V1, V2, V3 = (np.zeros((L, S), order='C'),
                      np.zeros((J, S), order='C'),
                      np.zeros((J, S), order='C'))
    else:
        V1, V2, V3 = init
    U1, U2, U3 = V1.copy(), V2.copy(), V3.copy()
    V1_prev, V2_prev, V3_prev = (np.empty_like(V1, order='C'),
                                 np.empty_like(V2, order='C'),
                                 np.empty_like(V3, order='C'))

    if constr == 'none':
        gamma = 1.0/(eta*(np.linalg.norm(np.dot(B.T, B)) + 1.0))
    else:
        gamma = 1.0/(eta*(np.linalg.norm(np.dot(B.T, B)) + 2.0))
    t = 1.
    # --

    # GAPS values
    gaps = list()
    primals = list()
    duals = list()

    for n in xrange(int(maxN)):
        t_prev = t
        V1_prev[:], V2_prev[:], V3_prev[:] = V1, V2, V3
        Gamma_aux[:] = Theta - eta*(np.dot(B.T, U1) + U2 + U3)

        # Data fit
        V1[:] = (1./(1. + gamma)) * (U1 + gamma*(np.dot(B, Gamma_aux) - Y))

        # Hard constraint (positivity and ...)
        grad = U2 + gamma*Gamma_aux
        if constr == 'pos':
            #V2[:] = arg - pos_proj(arg)                          # Pos
            V2[:] = neg_proj(grad)                                # (direct)
        elif constr == 'l2':
            V2[:] = grad - apply_by_row(ball_proj, grad, gamma)   # L2
        elif constr == 'posl2':
            V2[:] = grad - apply_by_row(pos_l2_proj, grad, gamma*RADIUS) # Pos+L2
        elif constr == 'posbox':
            V2[:] = grad - apply_by_row(box_proj, grad, gamma*UBOUND) # Pos+BOX
        #else V2 still remain 0

        # L1^2 norm
        grad = np.asfortranarray(U3 + gamma*Gamma_aux)
        V3 = apply_by_col(prox_l1_squared_norm, grad/gamma, tau/gamma)       
        V3[:] = grad - gamma * V3
        
        Gamma[:] = Theta - eta*(np.dot(B.T, V1) + V2 + V3)
        
        if constr == 'pos':
            PGamma[:] = apply_by_row(pos_proj, Gamma)           # Pos
        elif constr == 'l2':
            PGamma[:] = apply_by_row(ball_proj, Gamma)          # L2
        elif constr == 'posl2':
            PGamma[:] = apply_by_row(pos_l2_proj, Gamma, RADIUS)    # Pos+L2
        elif constr == 'posbox':
            PGamma[:] = apply_by_row(box_proj, Gamma, UBOUND)    # Pos+Box
        else:
            PGamma[:] = Gamma.copy()

        if n%10==0:
          primal = (
            (0.5 * np.sum((Y - np.dot(B, PGamma))**2)) +           # Data Fit
            (tau * (np.sum(np.sum(np.abs(PGamma), axis=0)**2))) +  # L1^2
            (np.sum((PGamma - Theta)**2) / (2.*eta))               # Prox
            # Projection not included because PGamma is feasible
            )
          dual = (
            ((np.sum(Gamma*Gamma) - ThetaNorm2) / (2.*eta)) +       # Prox*
            (0.5*np.sum(V1**2) + np.sum(V1*Y)) +                    # Fit*
            (np.sum(np.max(V3**2, axis=0))/(4.*tau))                # L1^2*
          )
          if constr == 'pos':
              pass
              # Pos: 0 because V2 is feasible (proj. onto the neg orthant)
          elif constr == 'l2':
              dual += np.sum(np.sqrt(np.sum(V2**2., axis=1)))         # ProjL2*
          elif constr == 'posl2':
              dual += RADIUS * np.sum(np.sqrt(np.sum(                 # ProjPos+L2*
                                    np.clip(V2, 0., np.inf)**2., axis=1)))
          elif constr == 'posbox':
              dual += UBOUND * np.sum(np.clip(V2, 0., np.inf))        # ProjPos+Box*
          #else constraint not included

          #print primal, dual

          gap = primal+dual
          gaps.append(gap)
          primals.append(primal)
          duals.append(dual)

        t = (1. + np.sqrt(1. + 4.*t_prev*t_prev)) * 0.5
        U1[:] = V1 + ((t_prev - 1.) / t) * (V1 - V1_prev)
        U2[:] = V2 + ((t_prev - 1.) / t) * (V2 - V2_prev)
        U3[:] = V3 + ((t_prev - 1.) / t) * (V3 - V3_prev)

        #if (n == 0) or ((n+1) % 1000 == 0):
        #    print 'Prox Phi #it: %d' % (n+1)
        #    print '   Gap: %.5f (th %.5e)' % (gap, (eps*eps)/(2.*eta))
        #    print '   Primal: %.5f' % primal
        #    print '   Dual: %.5f' % dual
        #    #t = 1. ####### Trick

        if gap <= (eps*eps)/(2.*eta):
            #print '##-Exit with gap %.5e (#%d it.)-##' % (gap, (n+1))
            break

    #if n == (maxN-1):
    #    print '##-Exit with gap %.5e (reached maximum #it %d)-##' % (gap, maxN)

    #print Theta
    #print Gamma
    #print PGamma
#
    return PGamma, gaps, primals, duals, (V1, V2, V3)

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
    mw = muw.ravel() # not diag...
    
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
        Dconj_U3[:] = apply_by_col(discr_derivate_conj, U3)
        Zeta_aux[:] = B - zeta*(np.dot(U1, Theta.T) + U2 + Dconj_U3)

        # Data fit
        V1[:] = (1./(1. + gamma)) * (U1 + gamma*(np.dot(Zeta_aux, Theta) - Y))

        # L1^2 norm
        grad = U2 + gamma*Zeta_aux
        V2 = apply_by_col(prox_l1_squared_norm, grad/gamma, lambda_/gamma)
        V2[:] = grad - gamma*V2
        
        # Weighted Total variation
        DZeta[:] = apply_by_col(discr_derivate, Zeta_aux)
        V3[:] = apply_by_col(interval_projection, U3 + gamma*DZeta, mw)

        Dconj_V3[:] = apply_by_col(discr_derivate_conj, V3)
        Zeta[:] = B - zeta*(np.dot(V1, Theta.T) + V2 + Dconj_V3)

        DZeta[:] = apply_by_col(discr_derivate, Zeta)

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

def cghDL(Y, J, lambda_, mu, tau, tvw=None, maxK=200, maxN=100, initB='pca',
          eps=1e-3):

    L, S = Y.shape
    init = initB

    maxiters = maxK
    precision = eps#1e-5

    #### B Initialization
    try:
        assert init.shape == (L, J)
        B0 = init.copy()
    except AttributeError:
        assert init in ['pca', 'rand']
        if init in 'pca':
            TMP = 1./np.sqrt(S-1) * Y.T
            TMP -= np.mean(TMP, axis=0)
            U, d, Vt = np.linalg.svd(TMP, full_matrices=False)
            V = np.array(Vt.T)
            B0 = V[:, :J]
        elif init in 'rand':
            sampling = np.arange(S)
            np.random.shuffle(sampling)
            B0 = Y[:,sampling[:J]]

    UP = 1.
    Theta0 = np.ones((J, S)) * UP

    PROJ = 'posbox'
    if tvw is None:
        w = np.ones((L-1, 1))        #### Total variation weigths
    else:
        w = tvw.reshape(L-1, 1)

    p = 2.                 #### Precision
    eta = 1./(S*J)
    zeta = 1./(L*J)

    tau *= (L/float(J))
    lambda_ *= (S/float(J))
    mu *= (S/float(J))

    #### Starting Duality Gap for PHI (with Vi=0, B fixed and Gamma=Theta0)
    PTheta0 = apply_by_row(box_proj, Theta0, UP)
    gap0_phi = (
        (0.5 * np.sum((Y - np.dot(B0, PTheta0))**2)) +
        (tau * (np.sum(np.sum(np.abs(PTheta0), axis=0)**2))) +
        (np.sum(PTheta0*PTheta0) / (2*eta)) +
        (np.sum(Theta0*Theta0) / (2*eta))
    )
    C_phi = np.sqrt(gap0_phi * 2.*eta)

    #### Starting Duality Gap for PSI (with Vi=0, Theta fixed and Zeta=B0)
    DB0 = np.apply_along_axis(discr_derivate, axis=0, arr=B0)
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

    B_diffs = list()
    Theta_diffs = list()
    B_prev = np.empty_like(B)
    Theta_prev = np.empty_like(Theta)

    for k in xrange(int(maxiters)):
        B_prev[:] = B.copy()
        Theta_prev[:] = Theta.copy()

        eps = 1. / ((k+1)**p)
        eta = 1. / ((k+1)**p)
        zeta = 1. / ((k+1)**p)
        
        (Theta, gaps,
         primals, duals,
         dual_var_phi) = prox_phi(Theta, eta, B, Y,
                                  tau, bound=UP,
                                  eps=C_phi*eps,
                                  maxN=maxN,
                                  init=dual_var_phi,
                                  constr=PROJ)#(k+1)**p)
        lastgapphi = gaps[len(gaps)-1]

        (B, gaps,
         primals, duals,
         dual_var_psi) = prox_psi(B, zeta, Theta, Y,
                                  mu*w, lambda_,
                                  eps=C_psi*eps,
                                  maxN=maxN,
                                  init=dual_var_psi)#(k+1)**p)
        
        lastgappsi = gaps[len(gaps)-1]

        B_diffs.append(np.sum((B - B_prev)**2)/np.sum(B_prev**2))
        Theta_diffs.append(np.sum((Theta - Theta_prev)**2)/np.sum(Theta_prev**2))

        convergence = (B_diffs[-1] <= precision and Theta_diffs[-1]<= precision)
        if convergence:
            break

    return {'B': B, 'Theta': Theta, 'conv': k,
            'gap_phi': lastgapphi, 'gap_psi': lastgappsi}

if __name__ == '__main__':
    np.random.seed(0)

    Y, ThetaTrue, BTrue = get_data(J=3, S=100, L=100)

    #### Parameters initialization
    lambda_ = 1e-1#mu**2         ### Atoms sparsity
    mu = 1e0                ### Atoms Total Variation
    tau = 1e-1             ### Weights sparsity

    print cghDL(Y, 10, lambda_, mu, tau ,tvw=None)['B'].sum()
