import numpy as np
import pylab as pl
from scipy import signal

from pycgh.datatypes.cytobands import CytoStructure, _chr2int
from pycgh.synth import ArrayCGHSynth
from pycgh.utils import array_trend
from pycgh.plots import profile, spatial
from pycgh.readers import agilent, ucsc_mapping

def global_median(acgh):
    chrX = acgh.F['chromosome'] == 23
    log2 = np.log2(acgh.F['test_signal']) - np.log2(acgh.F['reference_signal'])
    log2 -= np.median(log2[chrX]) # global median normalization
    return log2 # automatic defiltering

SYNTH = True
if SYNTH:
    cs = CytoStructure('data/ucsc/hg19/cytoBandIdeo.txt.gz')

    print 'Creating chip design'
    CHIP_DESIGN = ucsc_mapping('data/ucsc/hg19/agilentCgh4x44k.txt.gz',
                               filter_valid=True)

    print 'Creating data source'
    acgh_source = ArrayCGHSynth((430, 103), CHIP_DESIGN,
                                {'17q': [(4, 0.8), (2, 0.2)],
                                 '12': [(2, 0.8), (2, 0.2)],
                                 '2p': [(1, 0.8), (2, 0.2)],
                                 '1q': [(1, 0.8)]}, cs)
    print 'Created.'

    print 'Drawing synthetic CGH (median normalization, wrt X)'
    synth_acgh = acgh_source.draw()    
    synth_acgh['log2'] = global_median(synth_acgh)
    synth_acgh['true_log2'] = (np.log2(synth_acgh.F['true_test_signal']) -
                               np.log2(synth_acgh.F['true_reference_signal']))

    pl.figure()
    coords, cidx = profile(synth_acgh, signal=synth_acgh.F['log2'],
                           segmentation=synth_acgh.F['true_log2'])

    #pl.plot(coords, signal.medfilt(synth_acgh.F['log2'], 101), 'r-')
    #pl.plot(coords, synth_acgh.F['wave'], 'b-', alpha=0.5)
    pl.title('Synthetic')

    pl.figure()
    profile(synth_acgh, signal=synth_acgh.F['test_signal'],
            segmentation=signal.medfilt(synth_acgh.F['test_signal'], 101))
    a = pl.axis()

    pl.figure()
    profile(synth_acgh, signal=synth_acgh.F['reference_signal'],
            segmentation=signal.medfilt(synth_acgh.F['reference_signal'], 101))
    pl.axis(a)

    pl.figure()
    spatial(synth_acgh, signal=synth_acgh.F['log2'], median_center=True)
else:
    import os
    PATH = '/home/sabba/Phd/Tonini_IST/Dati/GEO/Samples/'
    #PATH = '/home/sabba/DISI/IST/GEO_Tonini/Samples/'

    print 'Reading Real CGH (median normalization wrt X)'
    PATH = os.path.join(PATH, '1562_G3.txt')
    #PATH = os.path.join(PATH, '1026_G2.txt')
    #PATH = os.path.join(PATH, '10003_G1_Cologne.txt')
    #PATH = os.path.join(PATH, '3383_G3_Cologne.txt')
    real_acgh = agilent(PATH, test_channel='r',
                        ucsc_mapping='data/ucsc/hg19/agilentCgh4x44k.txt.gz')
    real_acgh['log2'] = global_median(real_acgh)

    pl.figure()
    profile(real_acgh, signal=real_acgh.F['log2'],
            segmentation=signal.medfilt(real_acgh.F['log2'], 101))
    pl.title('Real')

    pl.figure()
    profile(real_acgh, signal=real_acgh.F['test_signal'],
            segmentation=signal.medfilt(real_acgh.F['test_signal'], 101))
    a = pl.axis()

    pl.figure()
    profile(real_acgh, signal=real_acgh.F['reference_signal'],
            segmentation=signal.medfilt(real_acgh.F['reference_signal'], 101))
    pl.axis(a)

    pl.figure()
    spatial(real_acgh, signal=real_acgh.F['log2'], median_center=True)

pl.show()
