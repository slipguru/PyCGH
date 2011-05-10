import sys

import numpy as np
from matplotlib import pylab as plt

import cghutils
from cghutils import plots

def main(file_name):
    # assumed test_channel='r' (as usual)
    aCGH = cghutils.AgilentCGH.load(file_name, fill_missings=True,
                                    qc_masking=True)

    #plt.figure()
    ratio = (np.log2(aCGH.filtered('test_signal')) -
             np.log2(aCGH.filtered('reference_signal')))
    #plots.spatial(ratio, aCGH.filtered('row'), aCGH.filtered('col'))

    #plt.figure()
    #plots.MA(aCGH.filtered('test_signal'), aCGH.filtered('reference_signal'))

    plt.figure()
    plots.profiles(aCGH)


if __name__ == '__main__':
    file_name = sys.argv[1]
    main(file_name)

    plt.show()
