import sys

import numpy as np
from matplotlib import pylab as plt

import cghutils
from cghutils import plots

def main(file_name):
    # assumed test_channel='r' (as usual)
    aCGH = cghutils.AgilentCGH.load(file_name, fill_missings=True,
                                    qc_masking=True)

    for fun in (plots.spatial, plots.MA, plots.profiles, plots.profile):
        plt.figure()
        fun(aCGH)

if __name__ == '__main__':
    file_name = sys.argv[1]
    main(file_name)

    plt.show()
