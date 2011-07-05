import os
from csv import DictReader

import numpy as np
from matplotlib import pylab as plt

from cghutils import AgilentCGH
from cghutils import plots as cghplots
from cghutils import r_functions as cghfunctions
#from cghutils.plots import array_image, MA_plot, cgh_profile
#from cghutils.r_functions import loess

ROOT_DIR = '/home/sabba/Phd/Tonini_IST'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'Group1')
SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')


# As example we pick the first GPL5477 sample reading the clinical data -------
clinical_info = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
for sample in clinical_info:
    if sample['Sample name'] == 'Sample 1':
        # Sample reading ======================================================
        # We have red=ch1=cy5 and green=ch2=cy3
        test_channel = 'r' if sample['ch1: source name'] != 'reference' else 'g'
        FILE_NAME = sample['Agilent Feature Extraction file']
        FILE_PATH = os.path.join(SAMPLES_DIR, FILE_NAME)
        aCGH = AgilentCGH.load(FILE_PATH, test_channel=test_channel,
                               qc_masking=True, fill_missings=False)


        # -- Ratios calculations ----------------------------------------------
        ratio = np.log2(aCGH['test_signal']) - np.log2(aCGH['reference_signal'])
        signal = 0.5 * (np.log2(aCGH['test_signal']) + np.log2(aCGH['reference_signal']))
        aCGH['ratio_median'] = (ratio-np.median(ratio))
        prediction = cghfunctions.loess(ratio, signal, span=0.03, degree=1,
                                        normalize=True,
                                        family='gaussian', iterations=3)
        aCGH['ratio_loess'] = (ratio-prediction)

        # -- M-A plot ---------------------------------------------------------
        plt.figure()
        plt.subplot(211)
        cghplots.MA(aCGH, title='MA plot - raw signals')
        plt.subplot(212)
        cghplots.MA(aCGH, M=(ratio-prediction), title='MA plot - loess')

        # -- Profiles ---------------------------------------------------------
        plt.figure()
        plt.subplot(211)
        cghplots.profile(aCGH, aCGH['ratio_median'])
        plt.subplot(212)
        cghplots.profile(aCGH, aCGH['ratio_loess'])


        plt.figure()
        pdf, bins, patches = plt.hist(aCGH['ratio_median'], 100, normed=True)
        plt.plot(0.5*(bins[1:] + bins[:-1]), pdf)


        # Sul mediano andrebbero messi dei paletti sui confini dei cromosomi

        from scipy import signal
        aCGH.sort(['chromosome', 'start_base'])
        plt.figure()
        plt.subplot(211)
        cghplots.profile(aCGH, signal.medfilt(aCGH['ratio_median'], 3))
        plt.subplot(212)
        cghplots.profile(aCGH, signal.medfilt(aCGH['ratio_loess'], 3))

        plt.figure()
        pdf, bins, patches = plt.hist(signal.medfilt(aCGH['ratio_median']),
                                      100, normed=True)
        plt.plot(0.5*(bins[1:] + bins[:-1]), pdf)

        #plt.figure()
        #pdf, bins = np.histogram(aCGH['ratio_median'], 100, normed=True)
        #print pdf
        #print bins
        #plt.plot(0.5*(bins[1:] + bins[:-1]), pdf)








        plt.show()

        #
        #fig.add_subplot(212)
        #MA_plot(acgh_valid['test_signal'], acgh_valid['ref_signal'],
        #        title='MA plot - raw signals agilent filtered')
        #
        ## -- Array Image: ref, test, ref_bg, ref_test -------------------------
        #fig = plt.figure()
        #fig.subplots_adjust(left=0.05, right=0.95)
        #loess_params = {'span': 1e-1, 'degree': 1, 'normalize': False,
        #                'family': 'symmetric', 'iterations': 4}
        #for i, k in enumerate(('test_signal', 'ref_signal', 'test_bg', 'ref_bg')):
        #    fig.add_subplot('32%d' % (i+1))
        #    out =  loess(acgh[k], acgh['row'], acgh['col'], **loess_params)
        #    array_image(out, acgh['row'], acgh['col'], title=k)
        #
        #fig.add_subplot('325')
        #out =  loess(acgh['test_signal'] - acgh['test_bg'],
        #             acgh['row'], acgh['col'], **loess_params)
        #array_image(out, acgh['row'], acgh['col'], title='test-bg')
        #fig.add_subplot('326')
        #out =  loess(acgh['ref_signal'] - acgh['ref_bg'],
        #             acgh['row'], acgh['col'], **loess_params)
        #array_image(out, acgh['row'], acgh['col'], title='ref-bg')
        #
        ## -- Array Image: raw ratio -------------------------------------------
        #fig = plt.figure()
        #raw_ratio = np.log2(acgh['test_signal']) - np.log2(acgh['ref_signal'])
        #out =  loess(raw_ratio, acgh['row'], acgh['col'], **loess_params)
        #fig.add_subplot('221')
        #array_image(out, acgh['row'], acgh['col'],
        #            title='raw log ratio (median centered)',
        #            median_center=True)
        #fig.add_subplot('222')
        #array_image(out[acgh['valid']], acgh_valid['row'], acgh_valid['col'],
        #            title='raw log ratio (median centered) - agilent filtered',
        #            median_center=True)
        #fig.add_subplot('223')
        #array_image(raw_ratio-out, acgh['row'], acgh['col'],
        #            title='raw log ratio (median centered)',
        #            median_center=True)
        #fig.add_subplot('224')
        #array_image((raw_ratio-out)[acgh['valid']], acgh_valid['row'], acgh_valid['col'],
        #            title='raw log ratio (median centered) - agilent filtered',
        #            median_center=True)


        exit()

