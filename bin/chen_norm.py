import os
from csv import DictReader

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab as ml

from mlabwrap import mlab
curr_dir = os.path.dirname(__file__)
mlab.addpath(os.path.join(curr_dir, '../libs/chen/'))
mlab.addpath(os.path.join(curr_dir, '../libs/chen/stprtool/'))
mlab.stprpath(os.path.join(curr_dir, '../libs/chen/stprtool/'))

from cghutils import ArrayCGH
from cghutils import plots as cghplots

ROOT_DIR = '/home/sabba/Phd/Tonini_IST'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'Group1')
SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')
OUT_DIR = os.path.join(ROOT_DIR, 'aCGH', 'Chen')

final_data = list()

# As example we pick the first GPL5477 sample reading the clinical data -------
clinical_info = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
for sample in clinical_info:
    #if sample['platform'] == 'GPL5477' and sample['Sample name'] != 'Sample 31':
    if sample['platform'] == 'GPL4093' and sample['Sample name'] != 'Sample 21':

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
                                                  0, nout=2)
        probes = mlab.cell2str(mlab.cellstr(aCGHdata.probe), '-').split('-') #Lento
        aCGH = ArrayCGH(probes,
                        np.asanyarray(aCGHdata.row.flatten(), dtype=int),
                        np.asanyarray(aCGHdata.col.flatten(), dtype=int),
                        aCGHdata.refInt.flatten(),
                        aCGHdata.sampleInt.flatten(),
                        np.asanyarray(aCGHdata.chrNumber.flatten(), dtype=int),
                        np.asanyarray(aCGHdata.start.flatten(), dtype=int),
                        np.asanyarray(aCGHdata.end.flatten(), dtype=int),
                        ratio_chen=aCGHdata.nRatioAdj.flatten())

        aCGH.save(os.path.join(OUT_DIR, sample['platform'], FILE_NAME))
