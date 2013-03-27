from itertools import chain

import numpy as np

# Rpy2 must be imported befor matplotlib to avoid a segmentation fault!
from utils import lowess, loess

from matplotlib import pylab as plt


def spatial(aCGH, signal=None, title=None, median_center=False,
            cmap=None, cbticks=None):
    plt.title('CGH Spatial plot' if title is None else title, size=8)

    rows = aCGH.F['row']
    cols = aCGH.F['col']

    if signal is None:
        signal = (np.log2(aCGH.F['test_signal']) -
                  np.log2(aCGH.F['reference_signal']))

    try:
        array_image = signal.reshape((rows.max(), cols.max()))
    except ValueError:
        array_image = np.ones((rows.max(), cols.max())) * np.nan
        for r, c, sig in zip(rows, cols, signal):
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


def MA(aCGH, M=None, title=None,
       points_color='k', median_color='b', lowess_color='r'):
    plt.title('CGH M-A plot' if title is None else title, size=8)

    test = aCGH.M['test_signal']
    reference = aCGH.M['reference_signal']
    A = (0.5 * (np.log2(test) + np.log2(reference))).compressed()

    if M is None:
        M = (np.log2(test) - np.log2(reference)).compressed()

    plt.scatter(A, M, label='Data', s=1, edgecolors='none', c=points_color)
    plt.xlabel('A')
    plt.ylabel('M')
    plt.axhline(0, lw=1, c='k', ls='--')
    plt.axhline(np.median(M), lw=1, c=median_color, label='Data Median')

    lowess_curve = lowess(A, M)
    x = np.asarray(lowess_curve[0])
    y = np.asarray(lowess_curve[1])

    plt.plot(x, y, '%s-' % lowess_color,
             lw=1, label='Data Lowess Curve') # Lowess Line

    plt.legend(loc='best', prop={'size': 8}, fancybox=True, shadow=False)


def profile(aCGH, signal=None,
            #chromosomes=None,
            ymin=None, ymax=None, cmin=-1, cmax=1, cmap=plt.cm.RdYlBu_r):
    """
        - signal: 'custom' cgh values
                  len(signal) == len(aCGH)
        - chromosomes: a list of chromosomes to plot from 1 to 24 (23=X, 24=Y)
        - segmentation: signal to plot against the cgh signal,
                        len(segmentation) == len(signal)
        - ymin, ymax: y axis limits
        - cmin, cmax: scatter plot color limits
    """

    chr_map = dict((str(c), c) for c in xrange(1, 25))
    chr_map['X'] = 23
    chr_map['Y'] = 24

    inverse_chr_map = dict((c, str(c)) for c in xrange(1, 23))
    inverse_chr_map[23] = 'X'
    inverse_chr_map[24] = 'Y'

    ## Chromosomes Fitering
    #if not chromosomes is None:
    #    chr_filtered = np.unique(chr_map[str(chr).strip().upper()]
    #                             for chr in chromosomes)
    #
    #    chridx = np.array(np.zeros(len(aCGH.F['chromosome'])), dtype=bool)
    #    for chr in chr_filtered:
    #        np.logical_or(chridx, aCGH.F['chromosome'] == chr, chridx)
    #
    #    #aCGH = aCGH.F[chridx]
    #    chromosomes = sorted(chr_filtered)
    #
    #    # Signal and segmentation filtering
    #    if not signal is None:
    #        signal = signal[chridx]
    #    if not segmentation is None:
    #        segmentation = segmentation[chridx]
    #else:

    # Signal Calculation
    if signal is None:
        test_signal = aCGH.M['test_signal']
        reference_signal = aCGH.M['reference_signal']
        signal = (np.log2(test_signal) - np.log2(reference_signal)).compressed()
    elif isinstance(signal, str) or isinstance(signal, unicode):
        signal = aCGH.F[signal]
    elif np.ma.getmask(signal) is np.ma.nomask:
        if len(signal) > aCGH.size: # Not filtered signal wrt current aCGH
            signal = signal[~aCGH['mask']] # Manual extraction
    else:
        signal = signal.compressed()

    # Plot-Coordinates Calculation
    from matplotlib import mlab as ml
    positions = aCGH.F[['chromosome', 'start_base']] # A compressed copy
    chromosomes = np.unique(positions['chromosome'])

    # Calculation of the max start_base for each Chromosome as shift
    coords = positions['start_base']
    summary = ml.rec_groupby(positions,
                             groupby=('chromosome',),
                             stats=(('start_base', np.max, 'shifts'),))
    shifts = np.zeros(len(summary)+1)
    shifts[1:] = np.cumsum(summary['shifts'])

    # Applying shifts
    for i, chr in enumerate(chromosomes):
         coords[positions['chromosome'] == chr] += shifts[i]

    # Chromosomes ticks and separators
    ticks = (shifts[:-1] + shifts[1:])/2.0
    separators = shifts[1:-1]

    # CGH Scatter plot
    plt.scatter(coords, signal, c=signal, cmap=cmap,
                vmin=cmin, vmax=cmax, s=1, edgecolors='none')

    # Axis limits
    if ymax is None:
        ymax = max(abs(min(np.nanmin(signal), -1.1)),
                   max(np.nanmax(signal), 1.1))
    if ymin is None:
        ymin = -ymax
    plt.axis([coords.min(), coords.max(), ymin, ymax])

    # Visual effects
    plt.colorbar(extend='both')
    plt.axhline(0.0, lw=1, color='gray', ls='--')
    plt.axhline(np.log2(3/2.), lw=1, color='orange', ls='--') # gain
    plt.axhline(1.0, lw=1, color='red', ls='--')
    plt.axhline(-1.0, lw=1, color='blue', ls='--')

    for sep in separators:
        plt.axvline(sep-1, lw=1, color='gray', ls='-')

    label_ticks = ['Chr %s' % inverse_chr_map[k] for k in chromosomes]
    plt.xticks(ticks, label_ticks, rotation=90, size=8)
    #plt.tick_params(axis='x', direction='out', length=3, colors='black',
    #                labelbottom='on')

    return coords

def profiles(aCGH, signal=None, vmin=-1, vmax=1, cmap=None, *args, **kwargs):
    coords = np.empty_like(signal)

    for i, chr in enumerate(chain(range(1, 23), ('X', 'Y'))):
        plt.subplot(6, 4, (i+1))
        c, ci = profile(aCGH, signal=signal, chromosome=chr, vmin=vmin, vmax=vmax,
                        cmap=cmap, *args, **kwargs)
        coords[ci] = c

    return coords
