import numpy as np
import pylab as pl
from scipy import signal

from pycgh.datatypes.cytobands import CytoStructure, _chr2int
from pycgh.synth import ArrayCGHSynth
from pycgh.utils import array_trend
from pycgh.plots import profile, spatial

cs = CytoStructure('data/ucsc/hg19/cytoBandIdeo.txt.gz')

data = np.loadtxt('data/ucsc/hg19/agilentCgh4x44k.txt.gz', dtype='S')
CHIP_DESIGN = {
    'P01': (1, 1, 200),
    'P02': (1, 300, 400),
    'P03': (1, 500, 600),
    'P04': (1, 700, 800),
    'P05': (2, 1, 200),
}

print 'Creating chip design'
CHIP_DESIGN = dict()
for _, chr, sb, eb, id in data:
    try:
        chr = _chr2int(chr[3:])
        CHIP_DESIGN[id] = chr, int(sb), int(eb)
    except:
        pass

print len(CHIP_DESIGN)

print 'Creating data source'
acgh_source = ArrayCGHSynth((430, 103), CHIP_DESIGN,
                            {'17': [(3, 0.8), (2, 0.2)],
                             '2p': [(1, 0.8), (2, 0.2)]}, cs)
print 'Created.'

print 'Starting Draw'
acgh = acgh_source.draw()
print 'drown'

log2 = np.log2(acgh.F['test_signal']) - np.log2(acgh.F['reference_signal'])
log2 -= np.median(log2) # global median normalization
acgh['log2'] = log2 # automatic defiltering


pl.figure()
profile(acgh, acgh.F['test_signal'])

pl.figure()
profile(acgh, acgh.F['reference_signal'])

pl.figure()
acgh.sort()
coords, cidx = profile(acgh, acgh.F['log2'],
                       segmentation=signal.medfilt(acgh.F['log2'], 101))

 # acgh already sorted
pl.plot(coords, np.log2(acgh.F['true_test_signal'] / 2.), '-')

pl.figure()
signal = (np.log2(acgh.F['test_signal']) -
          np.log2(acgh.F['reference_signal']))
signal = array_trend(signal, acgh.F['col'], acgh.F['row'])
spatial(acgh, signal=signal)#acgh.F['trend'])

pl.figure()
spatial(acgh, signal=acgh.F['trend'])

pl.figure()
spatial(acgh, median_center=True)

pl.show()
