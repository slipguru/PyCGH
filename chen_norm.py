import os
from csv import DictReader

from cghutils.plots import cgh_profile

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab as ml

from mlabwrap import mlab
mlab.addpath('libs/chen/')
mlab.addpath('libs/chen/stprtool/')
mlab.stprpath('libs/chen/stprtool/')

ROOT_DIR = '/home/sabba/Phd/Tonini_IST'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'Group1')
SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')

final_data = list()

# As example we pick the first GPL5477 sample reading the clinical data -------
clinical_info = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
for sample in clinical_info:
    if sample['platform'] == 'GPL5477':
        # Sample reading ======================================================
        # We have red=ch1=cy5 and green=ch2=cy3
        test_channel = 'r' if sample['ch1: source name'] != 'reference' else 'g'
        FILE_NAME = sample['Agilent Feature Extraction file']
        FILE_PATH = os.path.join(SAMPLES_DIR, FILE_NAME)

        print FILE_NAME,
        print '... ',

        #---------------------------
        data = mlab.ReadAgilentResult(FILE_PATH)
        aCGHdata, model = mlab.aCGHNormalization3(data, 1,
                                                  (0 if test_channel == 'r' else 1),
                                                  1, nout=2)

        print 'data conversion'

        keys = zip(mlab.cell2str(aCGHdata.probe, '-').split('-'),
                   aCGHdata.chrNumber.flatten(),
                   aCGHdata.start.flatten())
        data = dict(zip(keys, aCGHdata.nRatioAdj.flatten()))


        with open('CGHresult/%s' % FILE_NAME, 'w') as outfile:
            outfile.write('ProbeID\tChromosome\tPosition\tLogRatio\n')

            from operator import itemgetter
            keys = sorted(keys, key=itemgetter(1, 2))
            for k in keys:
                outfile.write('%s\t%d\t%d\t' % tuple(k))
                outfile.write('%.10f\n' % data[k])

