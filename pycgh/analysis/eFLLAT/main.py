# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Optimization main functions """

import matplotlib
matplotlib.use('Agg')

import os
import gc
import sys
import logging
import collections
import datetime as dt

import numpy as np
import pylab as pl

from bic import BIC_search
from plotting import plot_bics, pdf_plots, png_plots
   
### Imports current parameters ------------------------------------------------
from config import *
    
if __name__ =='__main__':
    ### Data ------------------------------------------------------------------
    data = np.load(sys.argv[1])
    Y = data['Y_0.1']
    #Y = data['Y_1.0']
    chrbrks = data['chrbreaks']
    if params['tvw'] == 'auto':
        params['tvw'] = chrbrks
    
    ### Output ----------------------------------------------------------------
    OUTPUT_DIR = '%s/Results/%s_%s' % (OUTPUT,
                                       ALGORITHM.upper(),
                                       dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    PART_DIR = os.path.join(OUTPUT_DIR, 'partials')
    os.makedirs(PART_DIR)
    
    ### Logging and saving ----------------------------------------------------
    logger = logging.getLogger('cghDL')
    logger.setLevel(logging.DEBUG)
    logfilename = '%s_%s.log' % (ALGORITHM, dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    lfile = logging.FileHandler(os.path.join(OUTPUT_DIR, logfilename))
    lfile.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    lfile.setFormatter(formatter)
    logger.addHandler(lfile)

    infomsg = '* Experiment with J={:2d}, lambda={:10.5f}, mu={:10.5f}, tau={:10.5f}'
    convmsg = ('  - iterat.: {:<5d}\n'
               '  - gap phi: {:<10.5g} | avg int. iters: {:<6.2f} | C: {:<10.5g}\n'
               '  - gap psi: {:<10.5g} | avg int. iters: {:<6.2f} | C: {:<10.5g}\n'
               '  - BIC:     {:<10.5f}')

    BICS = collections.defaultdict(list) # collected by callback!
    # CALL-BACK FUNCTION: called after each CGHDL fit #########################
    def callback(result, J, lambda_, mu, tau, BIC, eps_jumps):
        # Logging messages ----------------------------------------------------
        msg = infomsg.format(J, lambda_, mu, tau)
        print msg
        logger.info(msg)

        int_iters = result['internal_iters']

        msg = convmsg.format(result['conv']+1,
                             result['gap_phi'], np.mean(int_iters['phi'])+1, result['C_phi'],
                             result['gap_psi'], np.mean(int_iters['psi'])+1, result['C_psi'],
                             BIC)
        print msg
        logger.info('\n'+msg)
        ## --------------------------------------------------------------------

        # Current name (based on parameters)
        name = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J, lambda_, mu, tau)

        # Saving current result
        np.savez_compressed(os.path.join(PART_DIR, '%s.npz' % name), **result)

        # Updating BIC collection
        BICS[J].append(BIC)

        # PDF plot of the current solution
        pl.clf()
        filename = os.path.join(OUTPUT_DIR, name)
        png_plots(filename, data['YP'], result['Theta'], result['B'],
                  J, lambda_, mu, tau,
                  chrbrks, eps_jumps) # <-- main function parameters (closure)
        pl.clf()
        pl.close()
        gc.collect()
    ###########################################################################
    
    params['callback'] = callback
    result = BIC_search(Y, **params)
    
    final_msg = ('Best values: \n'
                 '  - J: %(J)d \n'
                 '  - lambda: %(lambda)f \n'
                 '  - mu: %(mu)f \n'
                 '  - tau: %(tau)f \n'
                 '  - BIC: %(BIC)f \n')
    logger.info(final_msg % result)
    
    # Plot BIC values for all the tuples of parameters
    plot_bics(os.path.join(OUTPUT_DIR, 'BICS.png'),
              params['J_range'], params['lambda_range'],
              params['mu_range'], params['tau_range'], BICS)
    
    # Plot, again, best solution in png
    filename = os.path.join(OUTPUT_DIR, 'Result')
    png_plots(filename, data['YP'], result['Theta'], result['B'],
              result['J'], result['lambda'], result['mu'], result['tau'],
              chrbrks, params['eps_jumps'])
    
    # Saving best result
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'Result.npz'),
                        Y=Y, BICS=BICS, **result)
    