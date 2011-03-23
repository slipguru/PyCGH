from itertools import chain

import numpy as np

# Rpy2 must be imported befor matplotlib to avoid a segmentation fault!
from rpy2 import robjects
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt

def array_image(rows, cols, signals, *args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    plt.title('CGH Spatial plot')

    # np.inf is masked ad showed as white
    array_image = np.ones((rows.max(), cols.max())) * np.inf
    for r, c, sig in zip(rows, cols, signals):
        array_image[r-1, c-1] = sig

    print rows.max()

    plt.imshow(array_image.T, cmap=plt.cm.jet,
                              aspect='equal',
                              interpolation='nearest',
                              norm=None)
    #plt.xticks(range(0, rows.max(), 10))
    #plt.grid(True)
    plt.xlabel('Rows')
    plt.ylabel('Columns')
    plt.colorbar(orientation='horizontal')

    return fig

def MA_plot(A, M, lowess=True, label='Input Signal'):
    plt.title('CGH M-A plot')

    plt.scatter(A, M, label=label, s=4, edgecolors='none', c='y')
    plt.xlabel(r'$A = \frac{log_2(r \cdot g)}{2}$')
    plt.ylabel(r'$M = log_2(\frac{r}{g})$')
    plt.axhline(np.median(M), lw=2, c='y', label='%s Median' % label)

    # Manual lowess normalization
    if lowess:
        rlowess = robjects.r['lowess'] # Lowess from R
        lowess_curve = rlowess(A, M, f=2./3, iter=3) #std params
        x = np.asarray(lowess_curve[0])
        y = np.asarray(lowess_curve[1])
        plt.plot(x, y, 'r-', lw=2, label='%s Lowess Curve' % label) # Lowess Line

        # M must be updated sorted
        assert np.allclose(x, np.sort(A))
        sorted_idxs = np.argsort(A) # as 'x' in lowess_curve
        M_norm = np.empty_like(M)
        M_norm[sorted_idxs] = M[sorted_idxs] - y

        plt.scatter(A, M_norm, label='Lowess-Normalized %s' % label,
                    s=4, edgecolors='none', c='c')
        plt.axhline(np.median(M_norm), lw=2, c='c',
                    label='Lowess-Normalized %s Median' % label)

        lowess_curve = rlowess(A, M_norm, f=2./3, iter=3) #std params
        x = np.asarray(lowess_curve[0])
        y = np.asarray(lowess_curve[1])
        plt.plot(x, y, 'b-', lw=2,
                 label='Lowess-Normalized %s Lowess Curve' % label)
        plt.legend()

        return M_norm

def cgh_profile(positions, signal, separators):
    plt.title('CGH profile')

    cmap = plt.get_cmap('jet')
    plt.scatter(positions, signal, c=signal, cmap=cmap,
                s=8, edgecolors='none')

    plt.axis([separators.min(), separators.max(),
              signal.min(), signal.max()])
    plt.colorbar()

    for sep in separators:
        plt.axvline(sep-1, lw=1, color='gray', ls='-')

    plt.axhline(0.0, lw=1, color='gray', ls='--')
    plt.axhline(1.0, lw=1, color='red', ls='--')
    plt.axhline(-1.0, lw=1, color='blue', ls='--')

    ticks = (separators[:-1] + separators[1:])/2.0
    label_ticks = ['Chr %s' % k for k in chain(range(1, 23), ('X', 'Y'))]
    plt.xticks(ticks, label_ticks, rotation=90)
    plt.tick_params(axis='x', direction='out', length=3, colors='black',
                    labelsize='small', labelbottom='on')
