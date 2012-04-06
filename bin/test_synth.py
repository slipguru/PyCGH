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
    print acgh.F['test_signal'].std()# / acgh.F['test_signal'].mean()
    print acgh.F['reference_signal'].std()# / acgh.F['reference_signal'].mean()

    return log2 # automatic defiltering

SYNTH = not True
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
                                noise = 0.1,
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

#exit()

#pl.figure()
#coords, cidx = profile(acgh, signal=acgh.F['log2'],
#                       segmentation=signal.medfilt(acgh.F['log2'], 101))
#a = pl.axis()
#pl.axis([a[0], a[1], -1.5, 1.5])
#
#pl.figure()
#pl.scatter(coords, acgh.F['test_signal'], marker='.', c='k')
#pl.plot(coords, signal.medfilt(acgh.F['test_signal'], 101))
#a = pl.axis()
##a = [a[0], a[1], 0, 1000]
##pl.axis(a)
#
#pl.figure()
#pl.scatter(coords, acgh.F['reference_signal'], marker='.', c='k')
#pl.plot(coords, signal.medfilt(acgh.F['reference_signal'], 101))
#pl.axis(a)

pl.figure()
count, bins, ignored = pl.hist(acgh.F['reference_signal'], 100, normed=True)
print 'Reference'
mu = np.mean(np.log(acgh.F['reference_signal']))
sigma = np.std(np.log(acgh.F['reference_signal']))

x = np.linspace(0.0, max(bins), 10000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi)))
pl.plot(x, pdf, linewidth=2, color='r')
pl.axis('tight')
pl.show()

pl.figure()
pl.hist(acgh.F['test_signal'], bins=100, align='mid', normed=True)

pl.figure()
pl.hist(acgh.F['log2'], bins=100, align='mid', normed=True)

#pl.figure()
#coords, cidx = profile(acgh, signal=acgh.F['log2'],
#                       segmentation=signal.medfilt(acgh.F['log2'], 101))
#pl.title('Real')

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
