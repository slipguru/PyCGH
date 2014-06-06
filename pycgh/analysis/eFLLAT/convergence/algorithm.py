# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Config-based main algorithm """

import numpy as np

import utils
import proxf

def minimize(Y, J, lambda_, mu, initB='pca', initTheta=None,
             maxK=200, maxN=100, eps=1e-3, callback=None):

    L, S = Y.shape

    #### B Initialization
    try:
        assert initB.shape == (L, J)
        B0 = initB.copy()
    except AttributeError:
        assert initB in ['pca', 'rand']
        B0 = utils.initB(Y, J, initB)

    #### Theta Initialization
    try:
        assert initTheta.shape == (J, S)
        Theta0 = initTheta.copy()
    except AttributeError:
        #Theta0 = np.ones((J, S)) * theta_bound
        Theta0 = np.zeros((J, S)) # for comparison with original FLLat

    #### Precision
    p = 1.         # How to justify such a choice?
    eta = 1./(S*J)
    zeta = 1./(L*J)

    # Parameters normalization
    lambda_ *= (S/float(J))  # l, j -- /S
    mu *= (S/float(J))       # l, j -- /S

    #### Starting Duality Gap for PHI (with Vi=0, B fixed and Gamma=Theta0)
    C_phi = proxf.c_phi(Y, Theta0, B0, eta)

    #### Starting Duality Gap for PSI (with Vi=0, Theta fixed and Zeta=B0)
    C_psi = proxf.c_psi(Y, Theta0, B0, lambda_, mu, zeta)

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
        
        # Option Const
        eta = 1.
        zeta = 1.
        
        # Option DESCR
        #eta = 1. / (k+1)
        #zeta = 1. / (k+1)
        
        # Option INCR
        #eta = (k+1)
        #zeta = (k+1)

        # Theta update
        Theta, phi_gaps, dual_var_phi = proxf.prox_phi(Theta, eta, B, Y,
                                                      eps=C_phi*epsk,
                                                      maxN=maxN,
                                                      init=dual_var_phi)

        # B update
        B, psi_gaps, dual_var_psi = proxf.prox_psi(B, zeta, Theta, Y, mu,
                                                  lambda_,
                                                  eps=C_psi*epsk,
                                                  maxN=maxN,
                                                  init=dual_var_psi)
        
        # Updating collections
        callback(phi_gaps, psi_gaps, k, J, lambda_, mu)

        B_diff = np.sum((B - B_prev)**2)/np.sum(B_prev**2)
        if k == 0: #first iteration
            Theta_diff = np.inf # -> more than eps
        else:
            Theta_diff = np.sum((Theta - Theta_prev)**2)/np.sum(Theta_prev**2)

        # convergence
        #if (B_diff <= eps and Theta_diff <= eps):
        print (B_diff <= eps and Theta_diff <= eps)
        #    break

    return {'B': B, 'Theta': Theta, 'conv': k,
            'gap_phi': phi_gaps[-1], 'gap_psi': psi_gaps[-1]}
