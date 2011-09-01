import os
from collections import defaultdict
import random

import numpy as np
from numpy.lib import arraysetops

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri
robjects.activate()

from matplotlib import pylab as plt
from matplotlib import mlab as ml

from cghutils import ArrayCGH, UCSC
from cghutils import plots as cghplots
from cghutils.utils import CytoBands, probes_average, LabeledMatrix

# R importing
importr('PMA')
pmd = robjects.r['PMD']
pmd_cv = robjects.r['PMD.cv']

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
PMD_DIR = os.path.join(ROOT_DIR, 'PMD')

lb = LabeledMatrix.load(os.path.join(PMD_DIR, 'cytobands_matrix.txt'))

# Data filtering
data = lb.asarray()
valid_columns = list()
for i, name in enumerate(sorted(lb.names)):
    if not np.any(np.isnan(data[:, i])):
        valid_columns.append(i)
data_filtered = data[:, valid_columns]
data_names = [sorted(lb.names)[i] for i in valid_columns]
data_sample_names = lb.samples

n, m = data_filtered.shape

K = 3
it_num = 5
V_out = np.empty((it_num, m, K))

for it in range(it_num):
    temp = range(n)
    random.shuffle(temp)
    temp = np.array(temp)
    idx, idx_rest =  temp[:n/2], temp[n/2:]
    print idx, idx_rest, len(idx) + len(idx_rest)

    ########################################
    out = pmd_cv(data_filtered[idx, :], type='ordered',
                 sumabsu=np.linspace(1, np.sqrt(len(idx)), 25),
                 niter=200, nfolds=5, center=True, chrom=robjects.NULL,
                 **{'lambda':robjects.NULL})
    #print out
    ## Ritorna sempre un valore alto del  parametro per minimizzare l'errore,
    ## che produce una selezione di quasi tutti i sample, il che si minimizza
    ## ad una sparse pca con sparsita' smooth sulle colonne.
    ## Questo potrebbe essere imputabile al fatto che i sample sono omogenei,
    ## tutti 4S
    best = out.rx2('bestsumabsu')
    print best


    # sumabsu: l1 on row
    # lambda: fl on col (NULL means detected by data)
    # K: number of factors?
    # center: data centered
    # chrom: indicates groups
    chrom = [b[:2] for b in data_names]
    out = pmd(data_filtered[idx_rest, :], type='ordered', sumabsu=best, #2,
              niter=100, K=K, center=True, chrom=robjects.NULL,
              **{'lambda':robjects.NULL})

    U = np.asarray(out.rx2('u'))
    V = np.asarray(out.rx2('v'))
    D = np.asarray(out.rx2('d'))
    v_init = np.asarray(out.rx2('v.init'))
    meanx = np.asarray(out.rx2('meanx'))
    lam = np.asarray(out.rx2('lambda'))

    V_out[it, :, :] = V[:, :]

    plt.figure()
    d, K = V.shape
    for k in range(K):
        #print
        nonzero = np.nonzero(V[:,k].squeeze())[0]
        nonzerou = np.nonzero(U[:,k].squeeze())[0]
        print ('factor %d => seleted var: %d, '
               'selected sample: %d') % (k, len(nonzero), len(nonzerou))
        print nonzerou
        plt.plot(np.arange(d)[nonzero],
                 V[nonzero,k], '.', label='factor %d' % k)

    plt.legend(loc='lower left')

    # Plot limits
    limits = dict()
    for i, chr in enumerate(chrom):
        limits[chr] = i

    for l in limits.values():
        plt.axvline(l, ls='-', lw=1, c='gray')
    plt.axhline(0, ls='-', lw=1, c='gray')

from cPickle import dump
with file('v_out.pkl', 'w') as f:
    dump(V_out, f)

V_freq = V_out.copy()
V_freq[V_out != 0.0] = 1

V_freq = V_freq.sum(axis=0)/float(it_num)
np.savetxt('freq.txt', V_freq, fmt='%.3f')

plt.figure()
for k in range(K):
    plt.plot(np.arange(m),
             V_freq[:,k], '.', label='factor %d' % k)

for l in limits.values():
    plt.axvline(l, ls='-', lw=1, c='gray')
plt.axhline(0, ls='-', lw=1, c='gray')

plt.show()
