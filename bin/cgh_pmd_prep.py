import os
from collections import defaultdict

import numpy as np
from numpy.lib import arraysetops

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt
from matplotlib import mlab as ml

from cghutils import ArrayCGH, UCSC
from cghutils import plots as cghplots
from cghutils.utils import CytoBands, probes_average, LabeledMatrix

PLATFORMS = {'GPL4093': UCSC['hg19']['agilentCgh2x105k'],
             'GPL5477': UCSC['hg19']['agilentCgh4x44k']}

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
IN_DIR = os.path.join(ROOT_DIR, 'Chen')
OUT_DIR = os.path.join(ROOT_DIR, 'PMD')

cb = CytoBands(release='hg19')

# File list calculation
files = list()
for plat in PLATFORMS:
    plat_dir = os.path.join(IN_DIR, plat)
    files.extend(list((plat, f)
                      for f in os.listdir(plat_dir) if f.endswith('txt')))

import itertools as it
expected_cytobands = dict()
for chr in range(1, 25):
    expected_cytobands.update(("%s%s" % (str(chr).zfill(2), b), np.nan)
                                for b in cb[chr].sub_bands())

# Data matrix
lb = LabeledMatrix(sorted(expected_cytobands))

for file_idx, (plat, filename) in enumerate(files):
    # I dati esportati da chen sono ordinati
    plat_dir = os.path.join(IN_DIR, plat)
    aCGH = ArrayCGH.load(os.path.join(IN_DIR, os.path.join(plat_dir, filename)))

    release = PLATFORMS[plat]
    print release, len(aCGH)

    # Get mappings
    from cghutils._agilent import _split_mapping as sm
    rel_data = np.genfromtxt(release, dtype=None, usecols=(4, 1, 2, 3),
                             names=('id', 'chr', 'sb', 'eb'))
    rel_data['sb'] += 1 # ucsc starts from zero
    rel_map = dict((id, sm('%s:%d-%d' % (c, s, e))) for id, c, s, e in rel_data)

    # Get useful data
    probe_id = aCGH['id']
    ratio = aCGH['ratio_chen']

    # Calculate the median of the probes replicated
    median_replicates = probes_average(probe_id, ratio, avg_function=np.median)

    # Calculate cytobands for probes
    cytobands = dict()
    for id in median_replicates.keys():
        try:
            chr, sb, eb, mask = rel_map[id]

            if not mask:
                bands = cb[chr][sb:eb] # bisogna gestire probe su piu' bande (poche)
                if len(bands) != 1: print bands
                cytobands[id] = "%s%s" % (str(chr).zfill(2), bands[0])
            else: # not sure mapped in the new release
                del median_replicates[id]
        except KeyError: #umapped in the new release
            del median_replicates[id]

    # Now we need a new list of ratio coherent with a list of cytobands
    ratio_filtered = list()
    cytobands_filtered = list()
    for i, id in enumerate(probe_id):
        if id in cytobands:
            ratio_filtered.append(ratio[i])
            cytobands_filtered.append(cytobands[id])

    # Calculate the mean ratio over the cytobands
    cytobands_mean = expected_cytobands.copy()
    cytobands_mean.update(probes_average(cytobands_filtered,
                                         ratio_filtered,
                                         avg_function=np.mean))

    cytobands, cb_ratios = zip(*sorted(cytobands_mean.items()))
    #plt.plot(cb_ratios)
    #plt.show()

    lb.append(filename, cb_ratios)

lb.save(os,path.join(OUT_DIR, 'cytobands_matrix.txt'))

#plt.show()
