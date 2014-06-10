# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Config-based main algorithm """

# Dynamic import --------------------------------------------------------------

from config import ALGORITHM

#import imp
#proxf = imp.load_module(ALGORITHM, *imp.find_module(ALGORITHM))

import importlib
proxf_module_name = 'pycgh.analysis.eFLLAT.%s' % ALGORITHM
proxf = importlib.import_module(proxf_module_name)

# -----------------------------------------------------------------------------

import numpy as np

import utils

def minimize(Y, J, lambda_, mu, tau, theta_bound=1.0, tvw=None,
             initB='pca', initTheta=None, maxK=200, maxN=100,
             eps=1e-3, eps_gap='auto', step=1.0):

    L, S = Y.shape

    #### B Initialization
    try:
        assert initB.shape == (L, J)
        B0 = initB.copy()
    except AttributeError:
        assert initB in ['pca', 'rand']
        B0 = utils.initB(Y, J, initB)

    #### Theta Initialization
    theta_bound = float(theta_bound)
    try:
        assert initTheta.shape == (J, S)
        Theta0 = initTheta.copy()
    except AttributeError:
        #Theta0 = np.ones((J, S)) * theta_bound
        Theta0 = np.zeros((J, S)) # for comparison with original FLLat

    #### Total variation weigths
    if tvw is None:
        w = np.ones((L-1, 1))
    else:
        w = tvw.reshape(L-1, 1)

    #### Precision
    p = 1.        # How to justify such a choice?
    eta = 1./(S*J)
    zeta = 1./(L*J)

    # Parameters normalization
    tau *= (L/float(J))      # s, j -- /J
    lambda_ *= (S/float(J))  # l, j -- /S
    mu *= (S/float(J))       # l, j -- /S

    #### Starting Duality Gap for PHI (with Vi=0, B fixed and Gamma=Theta0)
    C_phi = proxf.c_phi(Y, Theta0, B0, tau, eta, theta_bound)

    #### Starting Duality Gap for PSI (with Vi=0, Theta fixed and Zeta=B0)
    C_psi = proxf.c_psi(Y, Theta0, B0, lambda_, mu*w, zeta)

    dual_var_phi = None
    dual_var_psi = None
    B, Theta = B0, Theta0

    B_diff = None
    Theta_diff = None
    B_prev = np.empty_like(B)
    Theta_prev = np.empty_like(Theta)
    
    internal_iters = {'phi': [], 'psi': []}

    for k in xrange(int(maxK)):
        B_prev = B.copy()
        Theta_prev = Theta.copy()

        if eps_gap == 'auto':
            epsk = {'phi': C_phi/((k+1)**p), 'psi': C_psi/((k+1)**p)}
        else:
            epsk = {'phi': eps_gap, 'psi': eps_gap}
        
        if step == 'decr':
            eta = 1. / (k+1)
            zeta = 1. / (k+1)
        elif step == 'incr':
            eta = (k+1)
            zeta = (k+1)
        else:
            eta = step 
            zeta = step       

        # Theta update
        Theta, gap_phi, dual_var_phi, n_phi = proxf.prox_phi(Theta, eta, B, Y, tau,
                                                         bound=theta_bound,
                                                         eps=epsk['phi'],
                                                         maxN=maxN,
                                                         init=dual_var_phi)

        # B update
        B, gap_psi, dual_var_psi, n_psi = proxf.prox_psi(B, zeta, Theta, Y, mu*w,
                                                     lambda_,
                                                     eps=epsk['psi'],
                                                     maxN=maxN,
                                                     init=dual_var_psi)
        
        internal_iters['phi'].append(n_phi)
        internal_iters['psi'].append(n_psi)

        B_diff = np.sum((B - B_prev)**2)/np.sum(B_prev**2)
        if k == 0: #first iteration
            Theta_diff = np.inf # -> more than eps
        else:
            Theta_diff = np.sum((Theta - Theta_prev)**2)/np.sum(Theta_prev**2)

        # convergence
        if (B_diff <= eps and Theta_diff <= eps):
            break

    return {'B': B, 'Theta': Theta, 'conv': k,
            'gap_phi': gap_phi, 'gap_psi': gap_psi,
            'internal_iters': internal_iters,
            'C_phi': C_phi, 'C_psi': C_psi}

### Standard FLLat implementation called from R -------------------------------
if ALGORITHM == 'fllat_nowak':
    from fllat_nowak import minimize # overwrite!!
