import random

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
C = (1., -1., 1.)   # livelli segnale
J = (100, 500, 800) # inizio alterazione
K = (50, 100, 150)  # lunghezza alterazione
P = 5               # numero variabili latenti

bounds_err = 5.0  # probe
signal_err = 1e-1

#X = np.zeros((N, M))
rD = np.zeros((M, P))

## Data generation ------------------------------------------------------------
print 'Generazione dati...'
X = np.load('synth_cgh_data.npy')
#for i, x in enumerate(X):
#    x[:] = np.random.normal(scale=signal_err, size=M)
#
#    if i%2 == 0:
#        str, end = 0, 2
#    #elif i%3 == 1:
#        #str, end = 1, 2 # solo centrale
#    else:
#        str, end = 1, 3
#
#    for c, j, k in zip(C, J, K)[str:end]:
#
#        j += int(np.random.normal(scale=bounds_err))
#        k += int(np.random.normal(scale=bounds_err))
#
#        x[j:j+k] += c
#np.save('synth_cgh_data.npy', X)

print 'Generazione dizionario reale...'
for p, c, j, k in zip((1, 0, 2), C, J, K):
    rD[j:j+k,p] += c

### FLLat ----------------------------------------------------------------------
importr('FLLat')
##FLLat_PVE = robjects.r['FLLat.PVE']
FLLat_BIC = robjects.r['FLLat.BIC']
FLLat = robjects.r['FLLat']
#
### J selection ----------------------------------------------------------------
##J_seq = robjects.IntVector(range(1, N+1))
##result_pve = FLLat_PVE(X.T, J_seq)
##PVEs = np.asarray(result_pve.rx2('PVEs'))
##pl.plot(PVEs, 'o')
#
## Best model -----------------------------------------------------------------
print 'Calcolo modello...'
result_bic = FLLat_BIC(X.T, J=P)
lam1 = 1. * result_bic.rx2('lam1')[0]
lam2 = 1. * result_bic.rx2('lam2')[0]

## FLLat result --
D = X[random.sample(xrange(N), P), :].T
result_fflat = FLLat(X.T, J=P, B=D, lam1=lam1, lam2=lam2, **{'maxiter.T':1})
D = np.asarray(result_fflat.rx2('Beta'))
#A = np.asarray(result_fflat.rx2('Theta'))
#Xest = np.dot(D, A).T

#for a in A:
    #print np.linalg.norm(a, 2)

### Plot data estimations ------------------------------------------------------
#pl.figure(1)
#pl.title('FLLat')
#for p, (x, xest) in enumerate(zip(X, Xest)):
#    pl.subplot(6, 5, p+1)
#    pl.plot(x, 'b.')
#    pl.plot(xest, 'r-', alpha=0.7)
#
### Plot dictionary ------------------------------------------------------------
#pl.figure(2)
#pl.title('FLLat')
#for p, d in enumerate(D.T):
#    pl.subplot(4, 2, p+1)
#    pl.plot(d, '.')
#
### Plot coefficients ----------------------------------------------------------
#pl.figure(3)
#pl.title('FLLat')
#pl.imshow(A, aspect='auto', interpolation='nearest')
#pl.colorbar()
#
#print 'Minimo FLLat', np.linalg.norm(X - Xest, 2)

A = np.zeros((P, N))
#D = rD.copy()

## Plot dictionary ------------------------------------------------------------
pl.figure(2)
pl.title('Starting dictionary')
for p, d in enumerate(D.T):
    pl.subplot(4, 2, p+1)
    pl.plot(d, '.')

## FLLat + SPAMS results --
for i in xrange(500):

    # SPAMS --
    dmean = D.mean(axis=0)
    Dc = D - dmean
    A = mlab.test_FistaTree(X.T, Dc, A) # coefficients recalculation
    norms = np.apply_along_axis(np.linalg.norm, 1, A)
    norms = np.array([n or 1.0 for n in norms])
    A /= norms[:, np.newaxis] # limiting weights

    dsort = np.argsort(A.sum(axis=1)/N)
    D = D[:, dsort[::-1][[0, 1, 3, 2, 4]]]

    # FLLat --
    # proviamo ad iterare anche su A per limitare la norma... hack!
    result_fflat = FLLat(X.T, J=P, B=D, T=A, lam1=lam1, lam2=lam2,
                         **{'maxiter.T':0, 'maxiter.B':1, 'maxiter':100}) # coefficients updates ???
    D = np.asarray(result_fflat.rx2('Beta'))
    #A = np.asarray(result_fflat.rx2('Theta')) # A limitato

Xest = np.dot(D, A).T

print 'Minimo FLLat+SPAMS', np.linalg.norm(X - Xest, 2)
print A
for a in A:
    print np.linalg.norm(a, 2)

#exit()

## Plot data estimations ------------------------------------------------------
pl.figure(4)
pl.title('FLLat + SPAMS')
for p, (x, xest) in enumerate(zip(X, Xest)):
    pl.subplot(6, 5, p+1)
    pl.plot(x, 'b.')
    pl.plot(xest, 'r-', alpha=0.7)

## Plot dictionary ------------------------------------------------------------
pl.figure(5)
pl.title('FLLat + SPAMS')
for p, d in enumerate(D.T):
    pl.subplot(4, 2, p+1)
    pl.plot(d, '.')

## Plot coefficients ----------------------------------------------------------
pl.figure(6)
pl.title('FLLat + SPAMS')
pl.imshow(A, aspect='auto', interpolation='nearest')
pl.colorbar()

### Clustering atomi -----------------------------------------------------------
#print 'Clustering atomi...'
#linkage = ward(A)
#pl.figure()
#dendrogram(linkage)
#
### Clustering samples ----------------------------------------------------------
#print 'Clustering samples...'
#linkage = ward(A.T)
#pl.figure()
#dendrogram(linkage)


#pl.imshow(X, aspect='auto')
#pl.colorbar()
pl.show()
