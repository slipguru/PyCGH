import os
from collections import defaultdict

import numpy as np
from numpy.lib import arraysetops

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt
from matplotlib import mlab as ml

from cghutils import ArrayCGH
from cghutils import plots as cghplots
from cghutils.utils import CytoBands

#PLATFORM = 'GPL4093'
PLATFORM = 'GPL5477'

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
IN_DIR = os.path.join(ROOT_DIR, 'Chen', PLATFORM)

cb = CytoBands(release='hg18')

files = list(f for f in os.listdir(IN_DIR) if f.endswith('txt'))
filename = files[0]

aCGH = ArrayCGH.load(os.path.join(IN_DIR, filename))
autosomi = (aCGH['chromosome'] < 23)
aCGH = aCGH[autosomi]

bands_counter = defaultdict(int)
for id, chr, sb, eb in zip(aCGH['id'],
                           aCGH['chromosome'],
                           aCGH['start_base'],
                           aCGH['end_base']):
    try:
        bands_counter[str(chr) + cb[chr][sb:eb][0]] += 1
    except:
        pass

print min(bands_counter.values())
print max(bands_counter.values())
sidx = np.argsort(bands_counter.values())
print np.array(bands_counter.values())[sidx]
print np.array(bands_counter.keys())[sidx]
