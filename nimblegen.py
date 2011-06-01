import os
from csv import DictReader

import numpy as np
from matplotlib import pylab as plt

from cghutils import NimblegenCGH
from cghutils import plots as cghplots
from cghutils import r_functions as cghfunctions

ROOT_DIR = '/home/sabba/Phd/Zoppoli_DIMI'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'cgh_data')
#SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'def files/annotation_cgh.csv')

#signals = []
#signals_c = []
#annotation = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
#for i, sample in enumerate(annotation):
#    if sample['Tissue'].strip() != "male genomic reference" and sample['Ploidy'] == '3':
#        print sample['Filename']
#        filename = os.path.join(SAMPLES_DIR, sample['Filename'])
#        filename_c = os.path.join(SAMPLES_DIR, sample['Filename'][:-8] + '635.pair')
#        signals.append(np.loadtxt(filename, skiprows=2, usecols=(8,)))
#        signals_c.append(np.loadtxt(filename_c, skiprows=2, usecols=(8,)))
#    #if i == 3: break
#
#for sig, c in ((signals, 'r'), (signals_c, 'b')):
#    sig = np.asarray(sig)
#    sign_mean = sig.mean(axis=0)
#    sign_med = np.median(sig, axis=0)
#    sign_std = sig.std(axis=0)
#
#    coords = np.arange(len(sign_mean))
#    lowess_curve = cghfunctions.lowess(coords, sign_med, f=0.1)
#    x = np.asarray(lowess_curve[0])
#    y = np.asarray(lowess_curve[1])
#    #plt.errorbar(x, y, yerr=sign_std, fmt='r-', lw=2)
#    plt.plot(x, y, '%s-' % c, lw=2)
#
#plt.show()
#exit()

# Sample reading ======================================================
print 'Reading...'
#test_path = os.path.join(SAMPLES_DIR, '37968_532.pair')
#reference_path = os.path.join(SAMPLES_DIR, '37968_635.pair')
test_path = os.path.join(SAMPLES_DIR, '55650_532.pair')         # A549 tri
reference_path = os.path.join(SAMPLES_DIR, '55650_635.pair')
#test_path = os.path.join(SAMPLES_DIR, '59212_532.pair')         # HOP-62 htetr
#reference_path = os.path.join(SAMPLES_DIR, '59212_635.pair')
#test_path = os.path.join(SAMPLES_DIR, '59214_532.pair')         # HOP-92 ntetr
#reference_path = os.path.join(SAMPLES_DIR, '59214_635.pair')
#test_path = os.path.join(SAMPLES_DIR, '56192_532.pair')         # SF-295 pent
#reference_path = os.path.join(SAMPLES_DIR, '56192_635.pair')

#test_path = os.path.join(SAMPLES_DIR, '55668_635.pair')         # A549 tri
#reference_path = os.path.join(SAMPLES_DIR, '55668_532.pair')

aCGH = NimblegenCGH.load(test_path=test_path,
                         reference_path=reference_path)


# -- Ratios calculations ----------------------------------------------
print 'Global median...'

ratio = np.log2(aCGH['test_signal']) - np.log2(aCGH['reference_signal'])# * 3./2.)
aCGH['ratio_median'] = (ratio)#-np.median(ratio))

# -- M-A plot ---------------------------------------------------------
print 'MA plots...'
#plt.figure()
#plt.subplot(211)
#cghplots.MA(aCGH, title='MA plot - raw signals')
#plt.subplot(212)
#cghplots.MA(aCGH, M=aCGH['ratio_median'], title='MA plot - median')

## -- Array Image -----------------------------------------------------
print 'Array Image...'
#plt.figure()
#cghplots.spatial(aCGH['ratio_median'], aCGH['row'], aCGH['col'])

# -- Profiles ---------------------------------------------------------
print 'Profile...'
#plt.figure()
#coords = cghplots.profile(aCGH, aCGH['ratio_median'])
#coords = np.array(coords, dtype=float)

# -- Signals profiles -------------------------------------------------
#plt.figure()
lowess_curve = cghfunctions.lowess(coords, aCGH['test_signal'], f=0.1)
x = np.asarray(lowess_curve[0])
y = np.asarray(lowess_curve[1])
plt.plot(x, y, '%s-' % 'r',
         lw=2, label='Data Lowess Curve') # Lowess Line
plt.axhline(np.median(y), lw=1, c='k', ls='--')
print np.median(y)

lowess_curve = cghfunctions.lowess(coords, aCGH['reference_signal'], f=0.1)
x = np.asarray(lowess_curve[0])
y = np.asarray(lowess_curve[1])
plt.plot(x, y, '%s-' % 'b',
         lw=2, label='Data Lowess Curve') # Lowess Line
plt.axhline(np.median(y), lw=1, c='k', ls='--')
print np.median(y)

plt.show()

#path = 'NCI60_cell_line_panel_Genetics_Branch_I.R.Kirsch_1.xml'