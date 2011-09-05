import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

from mlabwrap import mlab
mlab.addpath('../spams_release/build/')
mlab.addpath('../spams_release/test_release/')

from matplotlib import pylab as pl

from scipy.cluster.hierarchy import ward, dendrogram


## Data parameters ------------------------------------------------------------
N = 30              # numero sample
M = 1000            # numero probe
S = 3               # numero segmenti
C = (1., -1., 2.)   # livelli segnale
J = (100, 500, 800) # inizio alterazione
K = (50, 20, 150)   # lunghezza alterazione

bounds_err = 5.0
signal_err = 1e-1

X = np.zeros((N, M))

## Data generation ------------------------------------------------------------
print 'Generazione dati...'
for i, x in enumerate(X):
    x[:] = np.random.normal(scale=signal_err, size=M)

    if i%2:
        str, end = 0, 2
    else:
        str, end = 1, 3

    for c, j, k in zip(C, J, K)[str:end]:

        j += int(np.random.normal(scale=bounds_err))
        k += int(np.random.normal(scale=bounds_err))

        x[j:j+k] += c

## FLLat ----------------------------------------------------------------------
importr('FLLat')
#FLLat_PVE = robjects.r['FLLat.PVE']
FLLat_BIC = robjects.r['FLLat.BIC']
FLLat = robjects.r['FLLat']

## J selection ----------------------------------------------------------------
#J_seq = robjects.IntVector(range(1, N+1))
#result_pve = FLLat_PVE(X.T, J_seq)
#PVEs = np.asarray(result_pve.rx2('PVEs'))
#pl.plot(PVEs, 'o')
J = 10

## Best model -----------------------------------------------------------------
print 'Calcolo modello...'
result_bic = FLLat_BIC(X.T, J=J, **{'maxiter.T':0}) #nessuna iterazione coeffs
lam1 = result_bic.rx2('lam1')
lam2 = result_bic.rx2('lam2')

result_fflat = FLLat(X.T, J=J, lam1=lam1, lam2=lam2)#, **{'maxiter.T':0})
D = np.asarray(result_fflat.rx2('Beta'))
A = np.asarray(result_fflat.rx2('Theta'))

pl.figure()
for p, d in enumerate(D.T):
    pl.subplot(5, 2, p+1)
    pl.plot(d, '.')

## Clustering -----------------------------------------------------------------
print 'Clustering matrice pesi...'
linkage = ward(A.T)
pl.figure()
dendrogram(linkage)

pl.figure()
pl.imshow(A, aspect='auto', interpolation='nearest')
pl.colorbar()


## SPAMS ----------------------------------------------------------------------
# Online learning for matrix factorization and sparse coding - FL
# Tree-structured sum of l2-norms
# mexFistaTree(Y,X,W0,tree,param):
#   - Y contiene i dati -> X
#   - X e' il dizionario -> D
#   - W i coefficienti -> A
#
# tree ... vedi test_FistaTree.m
# param.loss='square';
# param.regul='tree-l2';
# A = mexFistaTree(X,D,A0,tree,param)
# non c'e' l'implementazione completa del metodo...
# e nemmeno la sola fase separata di calcolo del dizionario dato A
# se uso fflat con maxiter.T = 0... dovrebbe tenermi sempre fisso Theta=A
# va modificata la funzione fflat in modo da accettare old.T oltre a old.B
param = mlab.struct('loss', 'square',
                    'regul', 'tree-l2')
tree = mlab.struct('own_variables', [0, 0, 3, 5, 6, 6, 8, 9],
                   'N_own_variables', [0, 3, 2, 1, 0, 2, 1, 1],
                   'eta_g', [1, 1, 1, 2, 2, 2, 2.5, 2.5],
                   'groups', [[0, 0, 0, 0, 0, 0, 0, 0],
                              [1, 0, 0, 0, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0, 0, 0, 0],
                              [1, 0, 0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 1, 0, 0, 0],
                              [0, 0, 0, 0, 1, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0, 1, 0]])

A, optim_info = mlab.mexFistaTree(X.T,D,A,tree,param)


#pl.imshow(X, aspect='auto')
#pl.colorbar()
pl.show()
