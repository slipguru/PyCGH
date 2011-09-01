import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

from matplotlib import pylab as pl


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
FLLat_PVE = robjects.r['FLLat.PVE']
FLLat_BIC = robjects.r['FLLat.BIC']

## J selection ----------------------------------------------------------------
#J_seq = robjects.IntVector(range(1, N+1))
#result_pve = FLLat_PVE(X.T, J_seq)
#PVEs = np.asarray(result_pve.rx2('PVEs'))
#pl.plot(PVEs, 'o')
J = 8

## Best model -----------------------------------------------------------------
result_bic = FLLat_BIC(X.T, J=J)
lam1 = result_bic.rx2('lam1')
lam2 = result_bic.rx2('lam2')
D = np.asarray(result_bic.rx2('opt.FLLat').rx2('Beta'))
A = np.asarray(result_bic.rx2('opt.FLLat').rx2('Theta'))

pl.figure()
for p, d in enumerate(D.T):
    pl.subplot(4, 2, p+1)
    pl.plot(d, '.')

pl.figure()
pl.imshow(A, aspect='auto', interpolation='nearest')
pl.colorbar()



#pl.imshow(X, aspect='auto')
#pl.colorbar()
pl.show()
