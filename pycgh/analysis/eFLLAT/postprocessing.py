# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" BIC-search related functions """

import itertools as it
import sys

import numpy as np

from bic import atoms_jumps
from plotting import png_plots, plot_breaks

import pylab as pl

from config import *

def dictionaries_plotting(J_range, lambda_range, mu_range, tau_range,
                          input_dir, output_dir):
    
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)

    # Parameters dimensions
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)

    for J in J_range:        
        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):
            # Current name (based on parameters)
            name = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J,
                                                                  lambda_range[j],
                                                                  mu_range[i],
                                                                  tau_range[k])
    
            # Loading saved result
            result = np.load(os.path.join(input_dir, '%s.npz' % name))

            B = result['B']
            Theta = result['Theta']
            
            # Dictionary plotting
            pl.figure()
            for idx, atom in enumerate(result['B'].T):
                pl.semilogy(np.abs(np.diff(atom)), '-', lw=2, label='#%d'%(idx+1))
            plot_breaks((np.where(chrbrks==0)[0]) if not chrbrks is None else [])
            pl.legend(loc='best')
            pl.savefig(os.path.join(output_dir, '%s_atoms_diff.png' % name))
            pl.clf()
            pl.close()
            
            pl.figure()
            for idx, atom in enumerate(result['B'].T):
                pl.semilogy(np.abs(atom), '-', lw=2, label='#%d'%(idx+1))
            plot_breaks((np.where(chrbrks==0)[0]) if not chrbrks is None else [])
            pl.legend(loc='best')
            pl.savefig(os.path.join(output_dir, '%s_atoms.png' % name))
            pl.clf()
            pl.close()

def offline_BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
                       input_dir, output_dir, eps_jumps=1e-3):
    L, S = Y.shape
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)

    # Parameters dimensions
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)
    
    # Initialization
    BIC_min = np.inf
    B_res, Theta_res = None, None
    J_res, mu_res, lambda_res, tau_res = None, None, None, None

    for J in J_range:        
        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):

            # Current name (based on parameters)
            name = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J,
                                                                  lambda_range[j],
                                                                  mu_range[i],
                                                                  tau_range[k])
    
            # Loading saved result
            result = np.load(os.path.join(input_dir, '%s.npz' % name))

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

    # Return the best result
    return {'B':B_res, 'Theta':Theta_res, 'J':J_res,
            'lambda':lambda_res, 'mu':mu_res, 'tau':tau_res, 'BIC':BIC_min}
    
    
def heuristic_BIC_search(Y, J_range, lambda_range, mu_range, tau_range,
                         input_dir, output_dir):
    L, S = Y.shape
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)

    # Parameters dimensions
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)
    
    # Initialization
    BIC_min = np.inf
    B_res, Theta_res = None, None
    J_res, mu_res, lambda_res, tau_res, eps_jumps_res = None, None, None, None, None

    for J in J_range:        
        # Evaluation
        for j, i, k in it.product(range(l), range(m), range(t)):

            # Current name (based on parameters)
            name = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J,
                                                                  lambda_range[j],
                                                                  mu_range[i],
                                                                  tau_range[k])
    
            # Loading saved result
            result = np.load(os.path.join(input_dir, '%s.npz' % name))

            B = result['B']
            Theta = result['Theta']
            
            # Autodetect eps_jumps
            auto_eps_dir = os.path.join(OUTPUT_DIR, 'Autodetected_EPS_Jumps', name)
            if not os.path.exists(auto_eps_dir):
                os.makedirs(auto_eps_dir)
                
            eps_jumps = autodetect_eps_jumps(B, auto_eps_dir)                       

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
                eps_jumps_res = eps_jumps

    # Return the best result
    return {'B':B_res, 'Theta':Theta_res, 'J':J_res,
            'lambda':lambda_res, 'mu':mu_res, 'tau':tau_res, 'BIC':BIC_min,
            'eps_jumps': eps_jumps_res}    
    
########## Autodetected Threshold
def epsmax(val):
    return max(np.finfo(float).eps, val)
    
def autodetect_eps_jumps(B, output_dir):
    eps_jumps = np.zeros(B.shape[1])
    
    for i, atom in enumerate(B.T):
        atomdiff = np.diff(atom)
        absatomdiff = np.abs(atomdiff)

        x = np.arange(50, 100, 0.1)
        percentiles = np.array([np.percentile(absatomdiff, q=xx) for xx in x])
        percentilesjumps = percentiles[1:]/percentiles[:-1]
        try:
            thidx = np.where(percentilesjumps >= 1e2)[0][0] # first one
        except IndexError:
            thidx = 0 # placeholder (median as threshold)
        
        th = percentiles[thidx]            
        eps_jumps[i] = th

        pl.figure()
        pl.semilogy(x, percentiles, '.-', lw=2, label=percentiles[thidx])
        pl.axvline(x[thidx], c='g')
        pl.title('Th = %.3e - %.2f percentile' % (th, x[thidx]))
        pl.savefig(os.path.join(output_dir, 'Atom_%d_TH.png' % (i+1)))
        
        # Filtering
        newatomdiff = atomdiff.copy()
        newatomdiff[absatomdiff <= th] = 0.0

        pl.figure()
        pl.subplot(2, 1, 1)
        pl.semilogy(absatomdiff, 'b-', label='Original ABS-Diff-Atom')
        pl.semilogy(map(epsmax, np.abs(newatomdiff)), 'r-', label='Filtered ABS-Diff-Atom')
        pl.axhline(th, c='g')
        pl.legend(loc='best', prop={'size':8})
        pl.xticks(np.arange(0, len(atomdiff), 10), rotation=45)
        pl.grid(True)

        pl.subplot(2, 1, 2)
        pl.plot(atom, 'b-', label='Original Atom')
        pl.plot(newatomdiff, 'r-', label='Filtered Diff-Atom')
        pl.legend(loc='best', prop={'size':8})
        pl.xticks(np.arange(0, len(atomdiff), 10), rotation=45)
        pl.grid(True)
    
        pl.savefig(os.path.join(output_dir, 'Atom_%d.png' % (i+1)))
        pl.clf()
        pl.close()
        
    return eps_jumps
    
    
def internal_iters(J_range, lambda_range, mu_range, tau_range, input_dir, output_dir):
    J_range = sorted(J_range)
    mu_range = sorted(mu_range)
    lambda_range = sorted(lambda_range)
    tau_range = sorted(tau_range)
    
    internal_iters_phi = []
    internal_iters_psi = []
    
    outfile = open(os.path.join(OUTPUT_DIR, 'saturations.txt'), 'w')

    for J in J_range:        
        # Evaluation
        for l, m, t in it.product(lambda_range, mu_range, tau_range):

            # Current name (based on parameters)
            name = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J, l, m, t)
    
            # Loading saved result
            result = np.load(os.path.join(input_dir, '%s.npz' % name))
            
            phi_iters = np.asarray(result['internal_iters'].tolist()['phi'], dtype=float)
            psi_iters = np.asarray(result['internal_iters'].tolist()['psi'], dtype=float)
            
            internal_iters_phi.append(phi_iters)
            internal_iters_psi.append(psi_iters)
            
            sat_phi = sat_psi = 'OK'
            outfile.write('%s\n' % name)
                       
            # Phi
            saturated = (phi_iters >= params['maxN']-1)            
            perc_saturated = (sum(saturated)/float(len(saturated))) * 100.;
            if np.any(saturated):
                sat_phi = '*** KO (%3.3f %%)***' % perc_saturated
            outfile.write('    - PHI: %s - ExitGap: %10.5e\n' % (sat_phi, result['gap_phi']))
            
            # Psi
            saturated = (psi_iters >= params['maxN']-1)
            perc_saturated = (sum(saturated)/float(len(saturated))) * 100.;
            if np.any(saturated):
                sat_psi = '*** KO (%3.3f %%)***' % perc_saturated
                
            outfile.write('    - PSI: %s - ExitGap: %10.5e\n' % (sat_psi, result['gap_psi']))
        
    outfile.close()  
    return {'phi': internal_iters_phi, 'psi': internal_iters_psi}
    
if __name__ == '__main__':
        
    PARTIALS_DIR = os.path.join(sys.argv[1], 'partials')
    OUTPUT_DIR = os.path.join(sys.argv[1], 'Results')
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    data = np.load(os.path.join(sys.argv[3]))
    Y = data['Y_0.1']
    chrbrks = data['chrbreaks']
    
    # Internal iterations
    internal_iters = internal_iters(params['J_range'], params['lambda_range'],
                                    params['mu_range'], params['tau_range'],
                                    PARTIALS_DIR, OUTPUT_DIR)
        
    pl.figure()
    pl.subplot(2, 1, 1)
    for i, step in enumerate(('phi', 'psi')):
        current_internal_iters = np.array(internal_iters[step]) + 1
        
        mean = np.mean(current_internal_iters, axis=0)
        median = np.median(current_internal_iters, axis=0)
        
        #mins = np.min(current_internal_iters, axis=0)
        #maxs = np.max(current_internal_iters, axis=0)
        mins = np.percentile(current_internal_iters, 10.0, axis=0)
        maxs = np.percentile(current_internal_iters, 90.0, axis=0)
                    
        pl.subplot(2, 1, i+1)
        pl.errorbar(np.arange(1, len(mean)+1), mean, yerr=np.c_[mean-mins, maxs-mean].T, label='mean')
        pl.plot(median, '-r', label='median')
        pl.title(step)
        pl.legend(loc='best')
        pl.axhline(params['maxN'], c='r', lw=2)
        pl.ylim(ymin=-params['maxN']*0.1, ymax=params['maxN']+(params['maxN']*0.1))
        pl.xlabel('External iteration (k)')
        pl.ylabel('Avg internal iterations (n)')
    
    pl.tight_layout()
    pl.savefig(os.path.join(OUTPUT_DIR, 'Iterations'))
    pl.clf()
    pl.close()
    # ----
    
    # Dictionaries
    dictionaries_plotting(params['J_range'], params['lambda_range'],
                          params['mu_range'], params['tau_range'],
                          PARTIALS_DIR, OUTPUT_DIR)
    # ----
    
    outfile = open(os.path.join(OUTPUT_DIR, 'BICS.txt'), 'w')
    msg = ('Best values: \n'
           '  - J: %(J)d \n'
           '  - lambda: %(lambda)f \n'
           '  - mu: %(mu)f \n'
           '  - tau: %(tau)f \n'
           '  - BIC: %(BIC)f \n')
    
    for eps_jumps in [0.0,
                      1e-17, 1e-16, 1e-15, 1e-14,
                      1e-13, 1e-12, 1e-11, 1e-10,
                      1e-9, 1e-8, 1e-7, 1e-6
                      ]:
        result = offline_BIC_search(Y, params['J_range'], params['lambda_range'],
                                    params['mu_range'], params['tau_range'],
                                    PARTIALS_DIR, OUTPUT_DIR, eps_jumps)
        
        header = '**** %.3e ****\n' % eps_jumps
        print header,
        print msg % result
        
        outfile.write(header)
        outfile.write(msg % result)
        
        # Plot, again, best solution in png
        filename = os.path.join(OUTPUT_DIR, 'Result_%.3e' % eps_jumps)
        png_plots(filename, data['YP'], result['Theta'], result['B'],
                  result['J'], result['lambda'], result['mu'], result['tau'],
                  chrbrks, eps_jumps)  
                
        # Saving best result
        np.savez_compressed(os.path.join(OUTPUT_DIR, 'Result_%.3e.npz' % eps_jumps),
                            Y=Y, **result)
    
    # Autodetected Threshold    
    result = heuristic_BIC_search(Y, params['J_range'], params['lambda_range'],
                                  params['mu_range'], params['tau_range'],
                                  PARTIALS_DIR, OUTPUT_DIR)
    eps_jumps = result['eps_jumps']
    header = '**** Autodetected ****\n'
    print header,
    print msg % result
    
    outfile.write(header)
    outfile.write(msg % result)
    
    # Plot, again, best solution in png
    filename = os.path.join(OUTPUT_DIR, 'Result_AUTO')
    png_plots(filename, data['YP'], result['Theta'], result['B'],
              result['J'], result['lambda'], result['mu'], result['tau'],
              chrbrks, eps_jumps)  
            
    # Saving best result
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'Result_AUTO.npz'),
                        Y=Y, **result)     
        
    outfile.close()