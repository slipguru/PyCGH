from itertools import chain

import numpy as np

# Rpy2 must be imported befor matplotlib to avoid a segmentation fault!
from r_functions import lowess, loess

from matplotlib import pylab as plt

def spatial(signals, rows, cols, title=None, median_center=False,
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


def MA(aCGH, M=None, title=None, points_color='k', median_color='b', lowess_color='r'):
    plt.title('CGH M-A plot' if title is None else title)

    test = aCGH['test_signal']
    reference = aCGH['reference_signal']
    nan_values = np.logical_or(np.isnan(test), np.isnan(reference))
    test = test[~nan_values]
    reference = reference[~nan_values]
    A = 0.5 * (np.log2(test) + np.log2(reference))

    if M is None:
        M = np.log2(test) - np.log2(reference)

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


def profile(aCGH, signal=None, chromosome=None, vmin=-1, vmax=1,
            cmap=None, title=None):
    #plt.title('CGH profile' if title is None else title)

    test_signal = aCGH['test_signal']
    reference_signal = aCGH['reference_signal']
    positions = aCGH[['chromosome', 'start_base']]

    if not chromosome is None:
        if chromosome == 'X': chr_val = 23
        elif chromosome == 'Y': chr_val = 24
        else: chr_val = int(chromosome)

        chridx = (aCGH['chromosome'] == chr_val)
        test_signal = test_signal[chridx]
        reference_signal = reference_signal[chridx]
        positions = positions[chridx]

    if signal is None:
        signal = np.log2(test_signal) - np.log2(reference_signal)

    # Calculation of the coordinates
    coords = positions['start_base']
    if chromosome is None:
        from matplotlib import mlab as ml
        summary = ml.rec_groupby(positions,
                                 groupby=('chromosome',),
                                 stats=(('start_base', np.max, 'shifts'),))
        shifts = np.array([0] + np.cumsum(summary['shifts']).tolist())


        for i in xrange(2, 25):
            coords[positions['chromosome'] == i] += shifts[i-1]

        ticks = (shifts[:-1] + shifts[1:])/2.0
        separators = shifts[1:-1]
    else:
        ticks = [coords.max() / 2.0]
    #--------------------------------

    cmap = plt.cm.jet if cmap is None else cmap
    plt.scatter(coords, signal, c=signal, cmap=cmap, vmin=vmin, vmax=vmax,
                s=8, edgecolors='none')

    max_h = max(abs(min(np.nanmin(signal), -1.1)), max(np.nanmax(signal), 1.1))
    plt.axis([coords.min(), coords.max(), -max_h, max_h])
    plt.colorbar()

    plt.axhline(0.0, lw=1, color='gray', ls='--')
    plt.axhline(1.0, lw=1, color='red', ls='--')
    plt.axhline(-1.0, lw=1, color='blue', ls='--')

    if chromosome is None:
        for sep in separators:
            plt.axvline(sep-1, lw=1, color='gray', ls='-')

        label_ticks = ['Chr %s' % k for k in chain(range(1, 23), ('X', 'Y'))]
        plt.xticks(ticks, label_ticks, rotation=90)
        plt.tick_params(axis='x', direction='out', length=3, colors='black',
                        labelsize='small', labelbottom='on')
    else:
        label_ticks = ['Chr %s' % chromosome]
        plt.xticks(ticks, label_ticks)
        plt.tick_params(axis='x', direction='out', length=3, colors='black',
                        labelsize='small', labelbottom='on')

    return coords

def profiles(aCGH, signal=None, vmin=-1, vmax=1, cmap=None):

    for i, chr in enumerate(chain(range(1, 23), ('X', 'Y'))):
        plt.subplot(6, 4, (i+1))
        profile(aCGH, signal=signal, chromosome=chr, vmin=vmin, vmax=vmax,
                cmap=cmap)
