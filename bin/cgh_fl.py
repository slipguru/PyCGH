import os

import numpy as np
from numpy.lib import arraysetops

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt
from matplotlib import mlab as ml

from cghutils import ArrayCGH
from cghutils import plots as cghplots

PLATFORM = 'GPL5477'

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
IN_DIR = os.path.join(ROOT_DIR, 'Chen', PLATFORM)
OUT_DIR = os.path.join(ROOT_DIR, 'FL', PLATFORM)

# R importing
importr('cghFLasso')
cghFLasso = robjects.r['cghFLasso']

for filename in list(f for f in os.listdir(IN_DIR) if f.endswith('txt'))[4:5]:

    # I dati esportati da chen sono ordinati
    aCGH = ArrayCGH.load(os.path.join(IN_DIR, filename))

    autosomi = (aCGH['chromosome'] < 23)
    out = cghFLasso(aCGH['ratio_chen'][autosomi],
                    chromosome=aCGH['chromosome'][autosomi],
                    nucleotide_position=aCGH['start_base'][autosomi],
                    FDR=1e-1)

    plt.figure()
    values = np.asanyarray(out.rx2('Esti.CopyN'))
    cghplots.profile(aCGH, signal=aCGH['ratio_chen'][autosomi],
                     indexes=autosomi,
                     # chromosome=17,
                     superimposed=values)

plt.show()
