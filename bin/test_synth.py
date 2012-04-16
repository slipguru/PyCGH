import numpy as np
import pylab as pl
from scipy import signal

from pycgh.datatypes.cytobands import CytoStructure, _chr2int
from pycgh.synth import ArrayCGHSynth
from pycgh.utils import array_trend, loess, lowess
from pycgh.plots import profile, spatial
from pycgh.readers import agilent, ucsc_mapping

def global_median(acgh):
    chrX = acgh.F['chromosome'] == 23
    log2 = np.log2(acgh.F['test_signal']) - np.log2(acgh.F['reference_signal'])
    log2 -= np.median(log2[chrX]) # global median normalization

    print log2[chrX].std()
    print log2[~chrX].std()

    return log2 # automatic defiltering

SYNTH = True
if SYNTH:
    cs = CytoStructure('data/ucsc/hg19/cytoBandIdeo.txt.gz')

    print 'Creating chip design'
    CHIP_DESIGN = ucsc_mapping('data/ucsc/hg19/agilentCgh4x44k.txt.gz',
                               filter_valid=True)

    print 'Creating data source'
    acgh_source = ArrayCGHSynth((430, 103), CHIP_DESIGN,
                                {
                                    '11q12': 3,
                                    '11q13': 3,
                                    '11q14': 1,
                                    '11q2': 1,
                                    '4': 1,
                                    '17q': 4
                                },
                                cytostructure = cs,
                                noise = 0.18, # doc: only an indication!
                                tissue_proportion = 0.5,
                                wave_bias_amplitude = 0.0)
    print 'Created.'

    print 'Drawing synthetic CGH (median normalization, wrt X)'
    acgh = acgh_source.draw('male')
    acgh['log2'] = global_median(acgh)

else:
    import os
    PATH = '/home/sabba/Phd/Tonini_IST/Dati/GEO/Samples/'
    #PATH = '/home/sabba/DISI/IST/GEO_Tonini/Samples/'

    print 'Reading Real CGH (median normalization wrt X)'
    #PATH = os.path.join(PATH, '1562_G3.txt')
    #PATH = os.path.join(PATH, '1026_G2.txt')
    #PATH = os.path.join(PATH, '10003_G1_Cologne.txt')
    #PATH = os.path.join(PATH, '3383_G3_Cologne.txt')
    PATH = os.path.join(PATH, '2216_G1.txt')
    acgh = agilent(PATH, test_channel='g',
                        ucsc_mapping='data/ucsc/hg19/agilentCgh4x44k.txt.gz')
    acgh['log2'] = global_median(acgh)


pl.figure()
pl.subplot(221)
coords = profile(acgh, signal=acgh.M['log2'])
pl.plot(coords, signal.medfilt(acgh.F['log2'], 101), c='gray')
pl.title('Log2 ratio')

pl.subplot(222)
pl.scatter(coords, acgh.F['test_signal'], marker='.', c='k')
pl.plot(coords, signal.medfilt(acgh.F['test_signal'], 101))
a = pl.axis()
a = (a[0], a[1], -300, 3700)
pl.axis(a)
pl.title('Raw Test Signal')

pl.subplot(223)
pl.scatter(coords, acgh.F['reference_signal'], marker='.', c='k')
pl.plot(coords, signal.medfilt(acgh.F['reference_signal'], 101))
pl.axis(a)
pl.title('Raw Reference Signal')

pl.subplot(224)
def lognormal(x, mu, sigma, alog=np.log):
    factor = alog(np.e)/(x*sigma*np.sqrt(2*np.pi))
    expo = np.exp(-0.5 * ( ((alog(x) - mu)*(alog(x) - mu)) / (sigma*sigma) ))
    return factor * expo

s = acgh.F['reference_signal']
logs = np.log2(s)
pl.hist(s, bins=1000, normed=True)
v = np.linspace(s.min(), s.max(), 10000)
pl.plot(v, lognormal(v, logs.mean(), logs.std(), np.log2), 'r-')
pl.title('Raw Reference Signal LogNormal Fit')

#### Normalizzazione cretina
#chrXY = acgh.F['chromosome'] >= 23
#r = acgh.F['reference_signal'][~chrXY]
#logr = np.log2(r)
#mean, std = logr.mean(), logr.std()
#ref_noise = r - 2**(mean + std*std*0.5) # "filtering" log-normal mean
#t = acgh.F['test_signal'][~chrXY] - ref_noise
#r = 2**(mean + std*std*0.5)
#
#log2 = np.log2(t) - np.log2(r)
#log2 -= np.median(log2) # global median normalization
#pl.figure()
#pl.scatter(coords[:len(t)], log2, marker='.', c=log2)
#pl.plot(coords[:len(t)], signal.medfilt(log2, 101))
#
#pl.figure()
#pl.scatter(coords[:len(t)], t, marker='.', c='k')
#pl.plot(coords[:len(t)], signal.medfilt(t, 101))
#pl.axis(a)






#################
#chr1 = (acgh.F['chromosome'] == 1)
#log2_loess = loess(acgh.F['log2'][chr1], coords[chr1], degree=2,
#                   span=0.1)

#log2_loess = loess(synth_acgh.F['log2'], coords, degree=1)
#pl.plot(coords, signal.medfilt(synth_acgh.F['log2'], 101), 'r-')
#pl.figure()
#pl.plot(coords[chr1], log2_loess, 'b-', lw=2)
#pl.scatter(coords[chr1], acgh.F['log2'][chr1], marker='.')

#pl.figure()
#profile(acgh, signal=acgh.F['test_signal'],
#        segmentation=signal.medfilt(acgh.F['test_signal'], 101))
#a = pl.axis()
#
#pl.figure()
#profile(acgh, signal=acgh.F['reference_signal'],
#        segmentation=signal.medfilt(acgh.F['reference_signal'], 101))
#pl.axis(a)
#
#pl.figure()
#spatial(acgh, signal=acgh.F['log2'], median_center=True)


#pl.figure()
#coords, cidx = profile(synth_acgh, signal=synth_acgh.F['log2'],
#                       segmentation=synth_acgh.F['true_log2'])
#
#chr1 = (synth_acgh.F['chromosome'] == 1)
#log2_loess = loess(synth_acgh.F['log2'][chr1], coords[chr1], degree=2,
#                   span=0.1)
#
##log2_loess = loess(synth_acgh.F['log2'], coords, degree=1)
##pl.plot(coords, signal.medfilt(synth_acgh.F['log2'], 101), 'r-')
#pl.figure()
#pl.plot(coords[chr1], log2_loess, 'b-', lw=2)
#pl.scatter(coords[chr1], synth_acgh.F['log2'][chr1], marker='.')
#
#pl.title('Synthetic')
#
#pl.figure()
#profile(synth_acgh, signal=synth_acgh.F['test_signal'],
#        segmentation=signal.medfilt(synth_acgh.F['test_signal'], 101))
#a = pl.axis()
#
#pl.figure()
#profile(synth_acgh, signal=synth_acgh.F['reference_signal'],
#        segmentation=signal.medfilt(synth_acgh.F['reference_signal'], 101))
#pl.axis(a)
#
#pl.figure()
#spatial(synth_acgh, signal=synth_acgh.F['log2'], median_center=True)

pl.show()
