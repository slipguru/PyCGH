import os
import cPickle as pkl

import numpy as np
from matplotlib import pylab as plt

from cghutils.utils import LabeledMatrix

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/RisultatiVari'
PMD_DIR = os.path.join(ROOT_DIR, 'PMD')


with file(os.path.join(PMD_DIR, 'v_out.pkl'), 'r') as f:
    V = pkl.load(f)

print V.shape

# Ogni colonna e' un fattore,
# le righe sono ordinate per citobande
# Il numero corrisponde alla frequenza

# Bisogna capire come interpretare i coefficienti prima
# di poter produrre qualcosa
