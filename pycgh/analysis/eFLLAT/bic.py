# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" BIC-search related functions """

import itertools as it

import numpy as np

import utils
import algorithm as alg 

def atoms_jumps(B, eps=1e-3):
    """ Counts the number of non-zero levels (up to eps tolerance). """
    ajumps = 0.0
    bs = np.empty_like(B[:,0])
    
    try:
        if len(eps) != B.shape[1]:
            raise ValueError('error in eps dimension, must be equal to the number of atoms')
    except TypeError:
        eps = np.ones(B.shape[1]) * eps

    for j in xrange(B.shape[1]):
        bs = B[:,j].copy()
        bs.sort()

        ajumps += (np.diff(bs) > eps[j]).sum() + 1.0

    return ajumps    
    
def BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
               theta_bound=1.0, tvw=None, initB='pca',
               maxK=200, maxN=100, eps=1e-3, eps_gap='auto', step=1.0,
               eps_jumps=1e-3, callback=None):

    L, S = Y.shape
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)

    # Parameters dimensions
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)

    # CGHDL params
    params = {'theta_bound':theta_bound, 'tvw':tvw,
              'maxK':maxK, 'maxN':maxN, 'eps':eps,
              'eps_gap':eps_gap, 'step': step}
    
    # Initialization
    BIC_min = np.inf
    B_res, Theta_res = None, None
    J_res, mu_res, lambda_res, tau_res = None, None, None, None

    for J in J_range:
        # For each new J value we start from following B0 and Theta0
        B0 = utils.initB(Y, J)
        #Theta0 = np.ones((J, S)) * theta_bound
        Theta0 = np.zeros((J, S)) # for comparison with original FLLat
        
        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):

            result = alg.minimize(Y, J,
                                  lambda_=lambda_range[j],
                                  mu=mu_range[i],
                                  tau=tau_range[k],
                                  initB=B0, initTheta=Theta0,
                                  **params)
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
                         lambda_range[j], mu_range[i], tau_range[k], BIC,
                         eps_jumps)

    # Return the best result
    return {'B':B_res, 'Theta':Theta_res, 'J':J_res,
            'lambda':lambda_res, 'mu':mu_res, 'tau':tau_res, 'BIC':BIC_min} 