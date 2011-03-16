import numpy as np
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


def MA_plot(red_signal, green_signal, *args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    plt.title('CGH M-A plot')

    A = 0.5 * (np.log2(red_signal) + np.log2(green_signal))
    M = np.log2(red_signal) - np.log2(green_signal)

    plt.scatter(A, M, s=4)
    plt.xlabel(r'$A = \frac{log_2(r \cdot g)}{2}$')
    plt.ylabel(r'$M = log_2(\frac{r}{g})$')
    plt.axhline(np.median(M), lw=2)

    return fig

def cgh_profile(signal, *args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    plt.title('CGH profile')

    plt.figure()
    plt.plot(signal, 'k.', lw=1)

    print xmin, xmax
    plt.axis([xmin, xmax, 0.5, 1.5])

    plt.axhline(0.0, lw=2)
    plt.axhline(1.0, lw=2)
    plt.axhline(-1.0, lw=2)


