import os

import numpy as np
from scipy.cluster.hierarchy import ward, dendrogram, leaves_list

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

from matplotlib import pylab as pl
from matplotlib import mlab as ml

import sys
sys.path.append('.')
from pycgh.datatypes.array_cgh import ArrayCGH
from pycgh.datatypes.cytobands import CytoBands
from pycgh.datatypes.dataset import Dataset
from pycgh.ucsc import UCSC

### FLLat ----------------------------------------------------------------------
importr('FLLat')
FLLat_PVE = robjects.r['FLLat.PVE']
FLLat_BIC = robjects.r['FLLat.BIC']
FLLat = robjects.r['FLLat']

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
FLLAT_DIR = os.path.join(ROOT_DIR, 'FLLat')

ds = Dataset.load(os.path.join(FLLAT_DIR, 'cytobands_matrix.txt'))

# Data filtering
data = ds.toarray()
valid_columns = list()
for i, name in enumerate(sorted(ds.names)):
    if not np.any(np.isnan(data[:, i])):
        valid_columns.append(i)
data_filtered = data[:, valid_columns]
data_names = [sorted(ds.names)[i] for i in valid_columns]
data_sample_names = ds.labels

chrom = [b[:2] for b in data_names]
limits = dict()
for i, chr in enumerate(chrom):
    limits[chr] = i #maximum index because are ordered

# Inputs ----------------------------------------------------------------------
X, (N, M) = data_filtered, data_filtered.shape

## J selection ----------------------------------------------------------------
#J_seq = robjects.IntVector(range(1, N+1))
#result_pve = FLLat_PVE(X.T, J_seq)
#PVEs = np.asarray(result_pve.rx2('PVEs'))
#pl.figure(2)
#pl.plot(PVEs, 'o')
P = 16

## Best model -----------------------------------------------------------------
print 'Calcolo parametri modello...',
result_bic = FLLat_BIC(X.T, J=P)
lam1 = result_bic.rx2('lam1')[0]
lam2 = result_bic.rx2('lam2')[0]
print lam1, lam2

## FLLat result ---------------------------------------------------------------
result_fflat = FLLat(X.T, J=P, B="pc", lam1=lam1, lam2=lam2)
D = np.asarray(result_fflat.rx2('Beta'))
A = np.asarray(result_fflat.rx2('Theta'))
Xest = np.dot(D, A).T

## Clustering atomi -----------------------------------------------------------
print 'Clustering atomi...'
linkage = ward(A)
pl.figure()
pl.title('Atoms clustering (by coefficients)')
dendrogram(linkage)
atoms_order = leaves_list(linkage)

## Clustering samples ---------------------------------------------------------
print 'Clustering samples...'
linkage = ward(A.T)
pl.figure()
pl.title('Samples clustering (by coefficients)')
dendrogram(linkage)
samples_order = leaves_list(linkage)

## Plots ----------------------------------------------------------------------

# Plot dictionary --
pl.figure()
pl.title('Dictionary Atoms (ordered by atoms clustering)')
vbound = min(np.abs(D.min()), np.abs(D.max()))
for p, d in enumerate(D[:,atoms_order].T):
    pl.subplot(4, 4, p+1)
    pl.title('Atom #%d' % atoms_order[p], size=10)
    pl.plot(d, '.')
    pl.axis([0, len(d), -vbound, vbound])

    ax = pl.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(10)

# Plot coefficients --
pl.figure()
pl.title('Coefficients Heatmap (ordered by clustering)')
vbound = min(np.abs(A.min()), np.abs(A.max()))
pl.imshow(A[atoms_order][:,samples_order],
          aspect='auto', interpolation='nearest', cmap=pl.cm.RdBu,
          vmin=-vbound, vmax=vbound)
pl.yticks(np.arange(len(atoms_order)), atoms_order)
pl.xticks(np.arange(len(samples_order)), samples_order)
pl.colorbar(extend='both')

# Plot raw data --
pl.figure()
pl.title('Raw data heatmap (ordered by samples clustering)')
vbound = min(np.abs(X.min()), np.abs(X.max()))
pl.imshow(X[samples_order],
          aspect='auto', interpolation='nearest', cmap=pl.cm.RdBu,
          vmin=-vbound, vmax=vbound)
for l in limits.values():
    pl.axvline(l, ls='-', lw=1, c='gray')
pl.yticks(np.arange(len(samples_order)), samples_order)
pl.colorbar(extend='both')

# Plot data estimations --
pl.figure()
pl.title('Estimated signals (ordered by samples clustering)')
vbound = min(np.abs(X.min()), np.abs(X.max()))
for p, (x, xest) in enumerate(zip(X[samples_order], Xest[samples_order])):
    pl.subplot(6, 6, p+1)
    pl.title('Sample #%d' % samples_order[p])
    pl.plot(x, 'b.')
    pl.plot(xest, 'r-', alpha=0.7)
    pl.axis([0, len(x), -vbound, vbound])

    ax = pl.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(10)




pl.show()
exit()














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
