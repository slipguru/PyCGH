from itertools import chain

import numpy as np

# Rpy2 must be imported befor matplotlib to avoid a segmentation fault!
from r_functions import lowess, loess

from matplotlib import pylab as plt

def array_image(signals, rows, cols, title=None, median_center=False,
                cmap=None, cbticks=None):
    plt.title('CGH Spatial plot' if title is None else title)

    try:
        array_image = signals.reshape((rows.max(), cols.max()))
    except ValueError:
        array_image = np.ones((rows.max(), cols.max())) * np.nan
        for r, c, sig in zip(rows, cols, signals):
            array_image[r-1, c-1] = sig

    if median_center:
        array_image -= np.median(array_image[~np.isnan(array_image)])
    vmin = np.nanmin(array_image)
    vmax = np.nanmax(array_image)

    xlabel, ylabel = 'columns', 'rows'
    r, c = array_image.shape
    if r > c:
        array_image = array_image.T
        xlabel, ylabel = ylabel, xlabel

    plt.tick_params(axis='both', direction='out', length=3, colors='black',
                    labelsize='small', labelbottom='on')
    plt.imshow(array_image,
               cmap=plt.cm.jet if cmap is None else cmap,
               aspect='equal',
               interpolation='nearest',
               norm=None,
               vmin=vmin, vmax=vmax)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    ticks = np.linspace(vmin, vmax, 10) if cbticks is None else np.asarray(cbticks)
    fmt = '%d' if ticks.dtype == np.int else '%.2f'
    cb = plt.colorbar(orientation='horizontal',
                      ticks=ticks, boundaries=cbticks,
                      aspect=40, format=fmt, spacing='proportional')
    cb.ax.tick_params(axis='x', direction='out', length=3, colors='black',
                      labelsize='small')

def MA_plot(test, reference, title=None,
            points_color='k', median_color='b', lowess_color='r'):
    plt.title('CGH M-A plot' if title is None else title)

    nan_values = np.logical_or(np.isnan(test), np.isnan(reference))
    test = test[~nan_values]
    reference = reference[~nan_values]

    M = np.log2(test) - np.log2(reference)
    A = 0.5 * (np.log2(test) + np.log2(reference))

    plt.scatter(A, M, label='Data', s=2, edgecolors='none', c=points_color)
    plt.xlabel('A')
    plt.ylabel('M')
    plt.axhline(0, lw=1, c=points_color, ls='--')
    plt.axhline(np.median(M), lw=2, c=median_color, label='Data Median')

    lowess_curve = lowess(A, M)
    x = np.asarray(lowess_curve[0])
    y = np.asarray(lowess_curve[1])

    plt.plot(x, y, '%s-' % lowess_color,
             lw=2, label='Data Lowess Curve') # Lowess Line

    plt.legend()


def cgh_profile(positions, signal, separators, vmin=-1, vmax=1):
    plt.title('CGH profile')

    cmap = plt.get_cmap('jet') #cool, jet
    plt.scatter(positions, signal, c=signal, cmap=cmap, vmin=vmin, vmax=vmax,
                s=8, edgecolors='none')

    plt.axis([positions.min(), positions.max(),
              np.nanmin(signal), np.nanmax(signal)])
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
