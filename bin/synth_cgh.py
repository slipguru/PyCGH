import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

from mlabwrap import mlab
mlab.addpath('../spams_release/build/')
mlab.addpath('../spams_release/test_release/')
mlab.addpath('bin/')

from matplotlib import pylab as pl

from scipy.cluster.hierarchy import ward, dendrogram


## Data parameters ------------------------------------------------------------
N = 30              # numero sample
M = 1000            # numero probe
S = 3               # numero segmenti
C = (1., -1., 2.)   # livelli segnale
J = (100, 500, 800) # inizio alterazione
K = (50, 200, 150)   # lunghezza alterazione

bounds_err = 5.0
signal_err = 1e-1

X = np.zeros((N, M))

## Data generation ------------------------------------------------------------
print 'Generazione dati...'
for i, x in enumerate(X):
    x[:] = np.random.normal(scale=signal_err, size=M)

    if i%2 == 0:
        str, end = 0, 2
    #elif i%3 == 1:
        #str, end = 1, 2 # solo centrale
    else:
        str, end = 1, 3

    for c, j, k in zip(C, J, K)[str:end]:

        j += int(np.random.normal(scale=bounds_err))
        k += int(np.random.normal(scale=bounds_err))

        x[j:j+k] += c


## Plot data ------------------------------------------------------------------
pl.figure()
for p, x in enumerate(X):
    pl.subplot(6, 5, p+1)
    pl.plot(x, '.')

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
P = 7

## Best model -----------------------------------------------------------------
print 'Calcolo modello...'
result_bic = FLLat_BIC(X.T, J=P, **{'maxiter.T':1}) #nessuna iterazione coeffs
lam1 = result_bic.rx2('lam1')
lam2 = result_bic.rx2('lam2')

print lam1, lam2

## Initial coefficients --
A = np.zeros((P, N))

## Fixed Dictionary (initial **guess**) --
#D = np.zeros((M, P))
#for c, j, k, i in zip(C, J, K, (1, 0, 4)):
    #D[j:j+k, i] += c
#D = X[np.random.random_integers(0, N-1, P), :].T # Random dictionary
result_fflat = FLLat(X.T, J=P, lam1=lam1, lam2=lam2, **{'maxiter.T':1})
D = np.asarray(result_fflat.rx2('Beta'))
A = np.asarray(result_fflat.rx2('Theta'))

### Algoritmo ------------------------------------------------------------------
#for i in xrange(100):
#
#    ## SPAMS --
#    A = mlab.test_FistaTree(X.T, D, A)
#
#    # FLLat --
#    result_fflat = FLLat(X.T, J=P, B=D, T=A, lam1=lam1, lam2=lam2, **{'maxiter.T':0})
#    D = np.asarray(result_fflat.rx2('Beta'))
#
#    #pl.figure()
#    #for p, d in enumerate(D.T):
#    #    pl.subplot(4, 2, p+1)
#    #    pl.plot(d, '.')

#A = mlab.test_FistaTree(X.T, D, A)
#print np.array(A > 1e-6, dtype=int)

## Plot dictionary ------------------------------------------------------------
pl.figure()
for p, d in enumerate(D.T):
    pl.subplot(4, 2, p+1)
    pl.plot(d, '.')

## Clustering atomi -----------------------------------------------------------
print 'Clustering atomi...'
linkage = ward(A)
pl.figure()
dendrogram(linkage)

## Clustering samples ----------------------------------------------------------
print 'Clustering samples...'
linkage = ward(A.T)
pl.figure()
dendrogram(linkage)

pl.figure()
pl.imshow(A, aspect='auto', interpolation='nearest')
pl.colorbar()



#pl.imshow(X, aspect='auto')
#pl.colorbar()
pl.show()
