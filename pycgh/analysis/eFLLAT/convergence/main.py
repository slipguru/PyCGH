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

import algorithm as alg

from plotting import pdf_plots
      
if __name__ =='__main__':
    ### Data ------------------------------------------------------------------
    data = np.load(sys.argv[1])
    Y = data['Y_0.1']
    
    ### Output ----------------------------------------------------------------
    OUTPUT_DIR = '%s/Results/Convergence_%s' % (os.path.expanduser('~'),
                                                dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    os.makedirs(OUTPUT_DIR)
    
    ### Logging and saving ----------------------------------------------------
    logger = logging.getLogger('cghDL')
    logger.setLevel(logging.DEBUG)
    logfilename = 'Convergence_%s.log' % dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    lfile = logging.FileHandler(os.path.join(OUTPUT_DIR, logfilename))
    lfile.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    lfile.setFormatter(formatter)
    logger.addHandler(lfile)

    infomsg = '* Experiment with J={:2d}, lambda={:10.5f}, mu={:10.5f}'
    convmsg = ('  - ext.it.: {:<5d}\n'
               '  - gap phi: {:<10.5g} | {:<10.5g}\n'
               '  - gap psi: {:<10.5g} | {:<10.5g}\n'
               '  - int.it.: phi->{:5d} | psi->{:5d}')

    PHI_GAPS = []
    PSI_GAPS = []
    ## CALL-BACK FUNCTION: called after each external iteration ################
    def callback(phi_gaps, psi_gaps, niter, J, lambda_, mu):
        # Logging messages ----------------------------------------------------
        msg = infomsg.format(J, lambda_, mu)
        print msg
        logger.info(msg)
    
        msg = convmsg.format(niter,
                             phi_gaps[0][-1], phi_gaps[1][-1],
                             psi_gaps[0][-1], psi_gaps[1][-1],
                             len(phi_gaps[0]), len(psi_gaps[0]))
        print msg
        logger.info('\n'+msg)
        ## --------------------------------------------------------------------
    
        # Updating collections
        PHI_GAPS.append(phi_gaps)
        PSI_GAPS.append(psi_gaps)
    ############################################################################
    
    J = 5
    lambda_ = 0.1
    mu = 0.1
    initB = 'pca'
    initTheta = None
    maxK = 10
    maxN = 1e4
    eps = -np.inf # not used!
    
    result = alg.minimize(Y, J, lambda_, mu, initB, initTheta,
                          maxK, maxN, eps, callback)
    
    B = result['B']
    Theta = result['Theta']
        
    ## Save results
    filename = os.path.join(OUTPUT_DIR, 'Result')
    np.savez_compressed('%s.npz' % filename,
                        Y=Y, B=B, Theta=Theta,
                        PHI_GAPS=PHI_GAPS,
                        PSI_GAPS=PSI_GAPS)
    
    fig = pl.figure(figsize=(10, 10))
    pl.suptitle(r'PSI_GAPS', weight='bold', size=8)
    for p, gaps in enumerate(PSI_GAPS):
        pl.subplot(5, 2, p+1)
        pl.plot(gaps[0], '.-b')
        pl.plot(gaps[1], '.-r')
    fig.tight_layout()
    pl.subplots_adjust(top=0.95)
    pl.savefig('%s_GAPS_PSI.png' % filename)
    
    fig = pl.figure(figsize=(10, 10))
    pl.suptitle(r'PSI_GAPS', weight='bold', size=8)
    for p, gaps in enumerate(PSI_GAPS):
        pl.subplot(5, 2, p+1)
        pl.semilogy(gaps[0], '.-b')
        pl.semilogy(gaps[1], '.-r')
    fig.tight_layout()
    pl.subplots_adjust(top=0.95)
    pl.savefig('%s_LogY_GAPS_PSI.png' % filename)
    
    fig = pl.figure(figsize=(10, 10))
    pl.suptitle(r'PHI_GAPS', weight='bold', size=8)
    for p, gaps in enumerate(PHI_GAPS):
        pl.subplot(5, 2, p+1)
        pl.plot(gaps[0], '.-b')
        pl.plot(gaps[1], '.-r')
    fig.tight_layout()
    pl.subplots_adjust(top=0.95)
    pl.savefig('%s_GAPS_PHI.png' % filename)
    
    fig = pl.figure(figsize=(10, 10))
    pl.suptitle(r'PHI_GAPS', weight='bold', size=8)
    for p, gaps in enumerate(PHI_GAPS):
        pl.subplot(5, 2, p+1)
        pl.semilogy(gaps[0], '.-b')
        pl.semilogy(gaps[1], '.-r')
    fig.tight_layout()
    pl.subplots_adjust(top=0.95)
    pl.savefig('%s_LogY_GAPS_PHI.png' % filename)
        
    # Plotting result
    pdf_plots('%s_result.pdf' % filename,
              data['Y'], Theta, B, J, lambda_, mu, 0.0, data['chrbreaks'], 1e-3)
              
    #pl.show()
        
    
    
    