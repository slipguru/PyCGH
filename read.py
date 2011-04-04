import os
from collections import defaultdict
from csv import DictReader

import numpy as np

from matplotlib import mlab
from cghutils.readers import AgilentReader, GPLReader
from cghutils.filters import probes_filter, split_mappings

ROOT_DIR = '/home/sabba/Phd/Tonini_IST'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'Group1')
SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')

# As example we pick the first GPL5477 sample reading the clinical data -------
clinical_info = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
for sample in clinical_info:
    #if sample['platform'] == 'GPL5477': break
    #if sample['Sample name'] == 'Sample 2': break

    # Saving csv for manor normalization
    if sample['Agilent Feature Extraction file'] == '1912_G1_251495014704_1_2.txt':
    #if sample['platform'] == 'GPL5477':
        # Sample reading ------------------------------------------------------
        # We assume red=ch1=cy5='reference' and green=ch2=cy3='tumor sample'
        swap = False if sample['ch1: source name'] == 'reference' else True
        FILE_NAME = sample['Agilent Feature Extraction file']
        FILE_PATH = os.path.join(SAMPLES_DIR, FILE_NAME)
        acgh = AgilentReader(FILE_PATH)

        # Remove controls and unmapped chr and extract data -------------------
        valid_probes = probes_filter(acgh)
        locations = acgh.feature('SystematicName')[valid_probes]
        probes = acgh.feature('ProbeName')[valid_probes]

        # Default: sample on the green=ch2=cy3 channel
        ch1_signal = acgh.feature('rMedianSignal')[valid_probes]
        ch2_signal = acgh.feature('gMedianSignal')[valid_probes]
        ch1_bg = acgh.feature('rBGMedianSignal')[valid_probes]
        ch2_bg = acgh.feature('gBGMedianSignal')[valid_probes]
        if sample['ch1: source name'] == 'reference':
            reference_signal, sample_signal = ch1_signal, ch2_signal
            reference_bg, sample_bg = ch1_bg, ch2_bg
        else:
            reference_signal, sample_signal = ch2_signal, ch1_signal
            reference_bg, sample_bg = ch2_bg, ch1_bg


        # Extact locus and sort by them -----------------------------------------------
        mappings = split_mappings(locations) # locus is a record array
        mappings_order = np.argsort(mappings) # sort by 'chromosome' and (if equal) by 'start'
                                            # (then, eventually, 'end')

        # Print Sample informations -------------------------------------------------
        print
        print '\n'.join(('Sample name: %(Sample name)s',
                         'Platform: %(platform)s')) % sample

        # Mappings fields are: 'chrmosome', 'start', 'end'
        summary = mlab.rec_groupby(mappings,
                                   groupby=('chromosome',),
                                   stats=(('chromosome', len, 'num_probes'),
                                          ('start', np.min, 'from'),
                                          ('end', np.max, 'to'),
                                          ))
        rjusts = (6, 10, 10, 10)
        print 'Chr # | # probes |   from   |    to    '
        print '------+----------+----------+----------'
        for record in summary:
            formatted = "|".join(str(val).rjust(just) for val, just in zip(record, rjusts))
            print formatted

        # M-A plot --------------------------------------------------------------------
        #plt.figure()
        A = 0.5 * (np.log2(sample_signal) + np.log2(reference_signal))
        M = np.log2(sample_signal) - np.log2(reference_signal)

        # Array Image plot ----------------------------------------------------------
        cols_max = acgh.feature('Col').max()
        rows_max = acgh.feature('Row').max()
        
        
        # Profile plot ----------------------------------------------------------------
        ratio = M

        # Version 1: only ordering of the positions
        positions = np.empty_like(ratio)
        separators = np.concatenate(([0], np.cumsum(summary['to'])))

        # Shifting...
        starting_locations = dict(zip(summary['chromosome'], separators[:-1]))
        for i, (chr, start, end) in enumerate(mappings):
            shifted_start = start + starting_locations[chr]
            shifted_end = end + starting_locations[chr]
            pos = int((shifted_start + shifted_end)/2) # median point of the probe

            positions[i] = pos

        # PROBLEMA: le probe replicate hanno stessa position...
        # lo scatter plot li stampa tutti (dorvrebbero essere visibili)


        break


