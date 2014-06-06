# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Data generation """

import random as rnd
import numpy as np
import pylab as pl
import datetime as dt
import os

from plotting import png_plots, plot_reconstruction

GAIN = np.log2(3./2.)
GAIN2 = np.log2(4./2.)
LOSS = np.log2(1./2.)

now = dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
OUTPUT_DIR = 'Datasets/%s' % now
os.makedirs(OUTPUT_DIR)

# 5 different events:
#   1. chr1 full half gain
#   2. chr2 small double gain
#   3. chr3 gain
#   4. chr3 loss
#   5. chr4 small half loss (border),
#      chr5 small half gain (border) + second small gain

# Each sample is a combination of 1, 2 or 3 atoms
# 1/5th of the sample are completely random

# Chromosomes have dimension 100: total 500 (=L)
# 5 samples for each event + 5 random: total 30 (=S) 

### Parameters ----------------------------------------------------------------
CHRNUM, CHRDIM = 5, 100
L, S, J = 500, 25, 5
B = np.empty((L, J))
indexes = np.array_split(np.arange(L), CHRNUM)

# TV breaks
chrbreaks = np.ones(L-1)
chrbreaks[CHRDIM-1::CHRDIM] = 0.0

### Atoms ---------------------------------------------------------------------
# Pattern 1
profile = np.zeros(L)
profile[indexes[0]] = GAIN / 2.
B[:, 0] = profile

# Pattern 2
profile = np.zeros(L)
profile[indexes[1][5:20]] = GAIN2
B[:, 1] = profile

# Pattern 3
profile = np.zeros(L)
profile[indexes[2]] = GAIN
B[:, 2] = profile

# Pattern 4
profile = np.zeros(L)
profile[indexes[2]] = LOSS
B[:, 3] = profile

# Pattern 5
profile = np.zeros(L)
profile[indexes[3][80:]] = LOSS / 2.
profile[indexes[4][:20]] = GAIN / 2.
profile[indexes[4][40:60]] = GAIN
B[:, 4] = profile

### Theta ---------------------------------------------------------------------
Theta = np.zeros((J, S))
atoms_indexes = range(J)
for s in xrange(S):
    
    # Avoid elision in order to get the simplest Theta solution
    rnd.shuffle(atoms_indexes)
    selected = atoms_indexes[:rnd.randint(1, 3)]
    while 2 in selected and 3 in selected:
        rnd.shuffle(atoms_indexes)
        selected = atoms_indexes[:rnd.randint(1, 3)]
                                 
    Theta[selected, s] = 1.0
    
### Data ----------------------------------------------------------------------
YT = np.dot(B, Theta)

# Perturbed Y
SP = S+5
YP = np.empty((L, SP))
YP[:, :S] = YT
YP[:, S:] = np.random.normal(loc=0.0, scale=0.5, size=(L, 5))

Y_0_01 = YP + np.random.normal(loc=0.0, scale=0.01, size=(L, SP))
Y_0_02 = YP + np.random.normal(loc=0.0, scale=0.02, size=(L, SP))
Y_0_05 = YP + np.random.normal(loc=0.0, scale=0.05, size=(L, SP))
Y_0_1  = YP + np.random.normal(loc=0.0, scale=0.1, size=(L, SP))
Y_0_2  = YP + np.random.normal(loc=0.0, scale=0.2, size=(L, SP))
Y_0_5  = YP + np.random.normal(loc=0.0, scale=0.5, size=(L, SP))
Y_1   = YP + np.random.normal(loc=0.0, scale=1.0, size=(L, SP))

png_plots('%s/Data' % OUTPUT_DIR, YT, Theta, B, J, 0.0, 0.0, 0.0, chrbreaks, 0.0)

breaks = np.where(chrbreaks==0)[0]
pl.figure()
plot_reconstruction(Y_0_01, YP, J, breaks, title='Noise 0.01')
pl.savefig('%s/Data_reconstruction_0.01.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_0_02, YP, J, breaks, title='Noise 0.02')
pl.savefig('%s/Data_reconstruction_0.02.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_0_05, YP, J, breaks, title='Noise 0.05')
pl.savefig('%s/Data_reconstruction_0.05.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_0_1, YP, J, breaks, title='Noise 0.1')
pl.savefig('%s/Data_reconstruction_0.1.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_0_2, YP, J, breaks, title='Noise 0.2')
pl.savefig('%s/Data_reconstruction_0.2.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_0_5, YP, J, breaks, title='Noise 0.5')
pl.savefig('%s/Data_reconstruction_0.5.png' % OUTPUT_DIR)
pl.figure()
plot_reconstruction(Y_1, YP, J, breaks, title='Noise 1.0')
pl.savefig('%s/Data_reconstruction_1.0.png' % OUTPUT_DIR)
pl.clf()
pl.close()

### Saving --------------------------------------------------------------------
np.savez('%s/dataset.npz' % OUTPUT_DIR,
                            **{
                                'Y': YT,
                                'YP': YP,
                                'B': B,
                                'Theta': Theta,
                                'Y_0.01': Y_0_01,
                                'Y_0.02': Y_0_02,
                                'Y_0.05': Y_0_05,
                                'Y_0.1':  Y_0_1,
                                'Y_0.2':  Y_0_2,
                                'Y_0.5':  Y_0_5,
                                'Y_1.0':  Y_1,
                                'chrbreaks':chrbreaks,
                                'CHRNUM': CHRNUM,
                                'CHRDIM': CHRDIM,
                                'CREATION': now
                            })
