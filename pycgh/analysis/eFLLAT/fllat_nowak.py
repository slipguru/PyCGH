# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" R FLLat wrapper """

import numpy as np

import utils

### R wrapping ----------------------------------------------------------------
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

importr('FLLat')
FLLat = robjects.r['FLLat']
### ---------------------------------------------------------------------------

## This function overwrite the standard "minimize" in algorithm
# The standard set of parameter is kept and ignored (eg. tau)
def minimize(Y, J, lambda_, mu, tau, theta_bound=1.0, tvw=None,
             initB='pca', initTheta=None, maxK=200, maxN=100, eps=1e-3,
             eps_gap='auto', step=1.0):

    L, S = Y.shape

    #### B Initialization
    try:
        assert initB.shape == (L, J)
        B0 = initB.copy()
    except AttributeError:
        assert initB in ['pca', 'rand']
        B0 = utils.initB(Y, J, initB)

    #### Theta Initialization
    # Theta0 = np.zeros((J, S)) # <- default

    # Parameters normalization
    #lambda_ *= (S/float(J))  # l, j -- /S
    #mu *= (S/float(J))       # l, j -- /S
    lambda_ *= 2.0*(S/float(J))  # l, j -- /S
    mu *= 2.0*(S/float(J))       # l, j -- /S
   
    result = FLLat(Y, J=J, B=B0, lam1=lambda_, lam2=mu,
                        thresh=eps,        # 1e-4 (default)
                        maxiter=maxK,      # 100  (default)
                        maxiter_B=maxN,    # 1    (  //   )  
                        maxiter_T=maxN)    # 1    (  //   )
    
    B = np.asarray(result.rx2('Beta'))
    Theta = np.asarray(result.rx2('Theta'))
    k = np.asarray(result.rx2('niter'))[0]
    gap_phi = -np.inf # placeholders
    gap_psi = -np.inf # --

    return {'B': B, 'Theta': Theta, 'conv': k,
            'gap_phi': gap_phi, 'gap_psi': gap_psi,
            'internal_iters': {'phi': [-1], 'psi': [-1]},
            'C_phi': -1, 'C_psi': -1}