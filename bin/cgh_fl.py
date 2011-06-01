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

PLATFORM = 'GPL4093'

ROOT_DIR = '/home/sabba/Phd/Tonini_IST/aCGH'
IN_DIR = os.path.join(ROOT_DIR, 'Chen', PLATFORM)
OUT_DIR = os.path.join(ROOT_DIR, 'FL', PLATFORM)

# R importing
importr('cghFLasso')
cghFLasso = robjects.r['cghFLasso']

TH = 0.2

cb = CytoBands(release='hg18')

gains_counter = defaultdict(int)
losses_counter = defaultdict(int)

files = list(f for f in os.listdir(IN_DIR) if f.endswith('txt'))
for filename in files:

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
    plt.savefig('out/%s.png' % filename)

    gains =  (values >= TH).flatten()
    losses = (values <= -TH).flatten()

    gains_set = set()
    for id, chr, sb, eb in zip(aCGH['id'][autosomi][gains],
                           aCGH['chromosome'][autosomi][gains],
                           aCGH['start_base'][autosomi][gains],
                           aCGH['end_base'][autosomi][gains]):
        try:
            band = cb[chr][sb:eb][0]
            #print chr, sb, eb, band, id
            #assert len(band) == 1
            gain = '%s%s' % (str(chr).zfill(2), band)
            gains_set.add(gain)
        except:
            pass

    for gain in gains_set:
        gains_counter[gain] += 1

    losses_set = set()
    for id, chr, sb, eb in zip(aCGH['id'][autosomi][losses],
                           aCGH['chromosome'][autosomi][losses],
                           aCGH['start_base'][autosomi][losses],
                           aCGH['end_base'][autosomi][losses]):
        try:
            band = cb[chr][sb:eb][0]
            #print chr, sb, eb, band, id
            #assert len(band) == 1

            loss = '%s%s' % (str(chr).zfill(2), band)
            losses_set.add(loss)
        except:
            pass

    for loss in losses_set:
        losses_counter[loss] += 1

import operator as op
gains_table = [(float(score)/len(files), name) for name, score in gains_counter.items()]
gains_table.sort(reverse=True)

losses_table = [(float(score)/len(files), name) for name, score in losses_counter.items()]
losses_table.sort(reverse=True)

with open('out/gains_%s.txt' % PLATFORM, 'w') as out:
    out.writelines('%s\t%s\n' % x for x in gains_table)

with open('out/losses_%s.txt' % PLATFORM, 'w') as out:
    out.writelines('%s\t%s\n' % x for x in losses_table)

#plt.show()
