# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" BIC-search related functions """

import itertools as it

import numpy as np

import utils
#import algorithm as alg 

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
    
#def BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
               #theta_bound=1.0, tvw=None, initB='pca',
               #maxK=200, maxN=100, eps=1e-3, eps_gap='auto', step=1.0,
               #eps_jumps=1e-3, callback=None):
               
def BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
               theta_bound=1.0, tvw=None,
               maxK=200, maxN=100, eps=1e-3, eps_gap='auto', step=1.0,
               eps_jumps=1e-3, callback=None):

    """
    This function runs the **E-FLLat** algorithm for all possible combination of parameters and selects the best model (i.e. tuple of parameters) using the Bayesian information criterion (BIC).
    
    Parameters
    ----------
    
    Y : numpy.ndarray
        The matrix representing all aCGH signals. Each column of the matrix represent a singla aCGH.
    
    J_range : array_like
        The list of all possible values for the parameter :math:`J`, representing the number of atoms which will be used for the reconstruction of the matrix.
        
    lambda_range : array_like
        All possible values for the :math:`\lambda` parameter.
    
    mu_range : array_like
        All possible values for the :math:`\mu` parameter.
    
    tau_range : array_like
        All possible values for the :math:`\tau` parameter.
    
    theta_bound : float, optional (default: ``1.0``)
        The value with which each entry of the :math:`\Theta` matrix will be initialized if it is not provided.
    
    tvw : array_like, optional (default: ``None``)
        An array containing the weights for the total variation term (when present).
    
    maxK : int, optional (default: 200)
        The maximum number of external iterations.
    
    maxN : int, optional (default: 100)
        The maximum number of internal iterations.
    
    eps : float, optional (default: 1e-3)
        The value used to stop the external iterations if neither of the two matrices changes more then this threshold.
    
    eps_gap : string, optional (default: ``"auto"``)
        The value used to stop the internal iterations when the solution does not change more than a fixed threshold between two iterations.
    
    step : float, optional (default: 1.0)
        A parameter controlling the step of internal iterations.
    
    eps_jumps : float, optional (default: 1e-3)
        The tolerance when computing the penalty for jumps in atoms.
    
    callback, function, optional (default: ``None``)
        If different from ``None``, the function is called at the end of each iteration and does something with the result obtained.
    """
    
    import algorithm as alg 

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