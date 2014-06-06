# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" BIC-search related functions """

import itertools as it

import numpy as np

import utils
import algorithm as alg 

    
def pseudo_BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
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
    

    for J in J_range:
        # For each new J value we start from following B0 and Theta0
        B0 = utils.initB(Y, J)
        #Theta0 = np.ones((J, S)) * theta_bound
        Theta0 = np.zeros((J, S)) # for comparison with original FLLat
        
        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):

            time = alg.minimize(Y, J,
                         lambda_=lambda_range[j],
                         mu=mu_range[i],
                         tau=tau_range[k],
                         initB=B0, initTheta=Theta0,
                         **params)
            
            # Callback execution
            if not callback is None:
                callback(time, J, lambda_range[j], mu_range[i], tau_range[k])
