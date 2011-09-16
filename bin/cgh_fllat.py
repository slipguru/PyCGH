import os
import csv

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

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH_GEO'
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
data_sample_names = [l[:-4] for l in ds.labels]

samples_colors = dict()
CSV = '/home/sabba/Phd/Tonini_IST/GEO/GSE255711_info.csv'
clinical_info = csv.DictReader(open(CSV, 'rb'), dialect=csv.excel_tab)
for sample in clinical_info:
    if sample['Sample_group'] == 'G1':
        c = 'r'
    elif sample['Sample_group'] == 'G2':
        c = 'b'
    elif sample['Sample_group'] == 'G3':
        c = 'g'
    else:
        assert False

    samples_colors[sample['Sample_title']] = c

# Che schifo!!!!
sample_labels = list()
sample_colors = list()
for name in data_sample_names:
    new_name = [n for n in samples_colors.keys() if name.startswith(n)]

    assert len(new_name) == 1
    sample_labels.append(new_name[0])
    sample_colors.append(samples_colors[new_name[0]])

chrom = [b[:2] for b in data_names]
limits = dict()
for i, chr in enumerate(chrom):
    limits[chr] = i #maximum index because are ordered

ordered_limits = sorted(limits.values())
chr_ticks = (ordered_limits + (np.array([0] + ordered_limits[:-1]))) / 2

# Inputs ----------------------------------------------------------------------
X, (N, M) = data_filtered, data_filtered.shape

## J selection ----------------------------------------------------------------
#J_seq = robjects.IntVector(range(1, N+1))
#result_pve = FLLat_PVE(X.T)#, J_seq)
#PVEs = np.asarray(result_pve.rx2('PVEs'))
#pl.figure(2)
#pl.plot(PVEs, 'o')
P = 16
#pl.show()
#exit()

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

## Plot function --------------------------------------------------------------

def plot_clustered_matrix(matrix, linkage_structure,
                          xticks, xlabels,
                          ylabels=None, ylabels_size=10,
                          ylimits=None,
                          left=0.2, width=0.7,
                          title=None,
                          colors=None,
                          fig=None):
    # Definitions for the axes
    if fig is None:
        fig = pl.figure()

    if ylabels is None:
        ylabels = 'Sample #%(num)d'

    bottom, height = 0.1, 0.8

    try:
        if len(linkage_structure) > 2:
            raise Exception('error')
        height -= 0.2
    except:
        linkage_structure = [linkage_structure]

    positions = (
        (0.05, bottom, left - 0.055, height), # left
        (left, bottom + height + 0.055, width, 0.2), # top
    )

    orientations = ('right', 'top')

    # Plot dendrograms (if one, only on left)
    orders = list()
    for ls, pos, ort in zip(linkage_structure, positions, orientations):
        axdend = fig.add_axes(pos)
        dendr = dendrogram(ls, orientation=ort)
        orders.append(dendr['leaves'])
        axdend.set_xticks([])
        axdend.set_yticks([])

    if len(orders) == 1:
        matrix = matrix[orders[0]]
    else:
        matrix = matrix[orders[0]][:, orders[1]]

    # Plot matrix
    axmat = fig.add_axes([left, bottom, width, height])
    vbound = min(np.abs(matrix.min()), np.abs(matrix.max()))
    im = axmat.matshow(matrix, aspect='auto',
                       cmap=pl.cm.RdBu_r,
                       vmin=-vbound, vmax=vbound)

    if ylimits:
        for l in ylimits.values():
            axmat.axvline(l, ls='-', lw=1, c='gray')
    axmat.set_xticks(xticks)
    axmat.set_xticklabels(xlabels,
                          size=10, rotation='vertical')
    axmat.set_yticks(np.arange(matrix.shape[0]))
    if title:
        pl.title(title, size=10, weight='bold')

    if type(ylabels) is str:
        axmat.set_yticklabels([ylabels % {'num':o} for o in orders[0]])
    else:
        axmat.set_yticklabels([ylabels[o] % {'num':o} for o in orders[0]])

    if colors:
        for tick, color in zip(axmat.yaxis.get_major_ticks(),
                              [colors[o] for o in orders[0]]):
            tick.label1On = False
            tick.label2On = True
            tick.label2.set_size(ylabels_size)
            tick.label2.set_color(color)
    else:
        for tick in axmat.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
            tick.label2.set_size(ylabels_size)


    # Plot colorbar
    axcolor = fig.add_axes([left, 0.05, width, 0.04])
    pl.colorbar(im, cax=axcolor, extend='both', orientation='horizontal')
    for tick in axcolor.get_xticklabels():
        tick.set_fontsize(10)

    return fig


## Clustering atomi -----------------------------------------------------------
print 'Clustering atomi...'
atoms_clustering = ward(A)
plot_clustered_matrix(D.T, atoms_clustering,
                      chr_ticks, ['Chr%d' % chr for chr in xrange(1, 25)],
                      ylabels='Atom #%(num)d', ylimits=limits,
                      title='Atoms ordered by coefficients clustering')

## Clustering samples ---------------------------------------------------------
print 'Clustering samples...'
#labels = ['%s [#%s]' % (name, '%(num)d') for name in data_sample_names]
samples_clustering = ward(A.T)
plot_clustered_matrix(X, samples_clustering,
                      chr_ticks, ['Chr%d' % chr for chr in xrange(1, 25)],
                      ylabels=sample_labels, ylimits=limits,
                      ylabels_size=7,
                      left=0.2,
                      width=0.6,
                      colors=sample_colors,
                      title='Samples ordered by coefficients clustering')

plot_clustered_matrix(Xest, samples_clustering,
                      chr_ticks, ['Chr%d' % chr for chr in xrange(1, 25)],
                      ylabels=sample_labels, ylimits=limits,
                      ylabels_size=7,
                      left=0.2,
                      width=0.6,
                      colors=sample_colors,
                      title='Samples estimation ordered by coefficients clustering')

# La generazione della figura va rivista!
samples_order = leaves_list(ward(A.T))
plot_clustered_matrix(A, [atoms_clustering, samples_clustering],
                      range(A.shape[1]), [str(a) for a in samples_order],
                      )

atoms_order = leaves_list(ward(A))
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


# Plot data estimations --
R = C = int(round(np.sqrt(N)))
if R*C < N:
    C += 1
pl.figure()
pl.title('Estimated signals (ordered by samples clustering)')
vbound = min(np.abs(X.min()), np.abs(X.max()))
for p, (x, xest) in enumerate(zip(X[samples_order], Xest[samples_order])):
    pl.subplot(R, C, p+1)
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
