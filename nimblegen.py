import os
from csv import DictReader

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt

from cghutils import NimblegenCGH
from cghutils import plots as cghplots
from cghutils import utils as cghfunctions

ROOT_DIR = '/home/sabba/Phd/Zoppoli_DIMI'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'cgh_data')
#SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'def files/annotation_cgh.csv')

# Sample reading ======================================================
print 'Reading...'
test_path = os.path.join(SAMPLES_DIR, '57428_532.pair')         # HCT-116 hdip
reference_path = os.path.join(SAMPLES_DIR, '57428_635.pair')
#test_path = os.path.join(SAMPLES_DIR, '55650_532.pair')         # A549 tri
#reference_path = os.path.join(SAMPLES_DIR, '55650_635.pair')
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

# -- Ploidy shifting ----------------------------------------------------------
#print 'Shifting...'
#ref_median = np.median(aCGH['reference_signal'])
#test_median = np.median(aCGH['test_signal'])
#aCGH['test_signal_shifted'] = aCGH['test_signal'] + (ref_median - test_median)
#aCGH['test_signal_shifted'] *= PLOIDY/2.

# So che la trisomia e' stata schiacciata al valore della mediana
#POLY = 4.
#mediana_attesa = POLY/2.
#segnale = aCGH['test_signal'] / aCGH['reference_signal']
#segnale_atteso = (mediana_attesa * segnale) / np.median(segnale)
#print np.median(segnale)

# -- Profiles ---------------------------------------------------------
print 'Profile...'
plt.figure()
coords = cghplots.profile(aCGH)#, signal=ratio)#(ratio-np.median(ratio))+center)
coords = np.array(coords, dtype=float)
#plt.axhline(center, lw=1, c='g', ls='--')

plt.axhline(np.log2(3./2.), lw=1, c='k', ls='-')
plt.axhline(np.log2(4./2.), lw=1, c='k', ls='-')
plt.axhline(np.log2(5./2.), lw=1, c='k', ls='-')
plt.axhline(np.log2(1./2.), lw=1, c='k', ls='-')


# -- Signals profiles -------------------------------------------------
#plt.figure()
#for signal in ('test_signal', 'reference_signal'):
#    lowess_curve = cghfunctions.lowess(coords, aCGH[signal], f=0.1)
#    x = np.asarray(lowess_curve[0])
#    y = np.asarray(lowess_curve[1])
#    plt.plot(x, y, '-',
#             lw=2, label=signal) # Lowess Line
#    plt.axhline(np.median(y), lw=1, c='k', ls='--')
#plt.legend()

# Shifted Profile
#plt.figure()
#ratio = np.log2(aCGH['test_signal_shifted']) - np.log2(aCGH['reference_signal'])
#aCGH['ratio_shifted'] = (ratio)#-np.median(ratio)
#cghplots.profile(aCGH, signal=aCGH['ratio_shifted'])
#plt.axhline(np.log2(PLOIDY/2.), lw=1, c='k', ls='--')
#plt.axhline(np.median(ratio), lw=1, c='g', ls='--')


# -- Ratios calculations ----------------------------------------------
#print 'Global median...'
#ratio = np.log2(aCGH['test_signal']) - np.log2(aCGH['reference_signal'])# * 3./2.)
#aCGH['ratio_median'] = (ratio)-np.median(ratio)

# -- M-A plot ---------------------------------------------------------
#print 'MA plots...'
#plt.figure()
#plt.subplot(211)
#cghplots.MA(aCGH, title='MA plot - raw signals')
#plt.subplot(212)
#cghplots.MA(aCGH, M=aCGH['ratio_median'], title='MA plot - median')

## -- Array Image -----------------------------------------------------
#print 'Array Image...'
#plt.figure()
#cghplots.spatial(aCGH['ratio_median'], aCGH['row'], aCGH['col'])



# -- Signals profiles -------------------------------------------------
#plt.figure()
#lowess_curve = cghfunctions.lowess(coords, aCGH['test_signal'], f=0.1)
#x = np.asarray(lowess_curve[0])
#y = np.asarray(lowess_curve[1])
#plt.plot(x, y, '%s-' % 'r',
#         lw=2, label='Data Lowess Curve') # Lowess Line
#plt.axhline(np.median(y), lw=1, c='k', ls='--')
#print np.median(y)
#
#lowess_curve = cghfunctions.lowess(coords, aCGH['reference_signal'], f=0.1)
#x = np.asarray(lowess_curve[0])
#y = np.asarray(lowess_curve[1])
#plt.plot(x, y, '%s-' % 'b',
#         lw=2, label='Data Lowess Curve') # Lowess Line
#plt.axhline(np.median(y), lw=1, c='k', ls='--')
#print np.median(y)

plt.show()

#path = 'NCI60_cell_line_panel_Genetics_Branch_I.R.Kirsch_1.xml'
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
