import os
from collections import defaultdict
from csv import DictReader

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab

from cghutils.readers import AgilentReader, GPLReader
from cghutils.filters import probes_filter, split_mappings
from cghutils.plots import array_image, MA_plot, cgh_profile

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
        M_norm = MA_plot(A, M, lowess=True, label='Raw Signal')

        # Array Image plot ----------------------------------------------------------
        cols = acgh.feature('Col')[valid_probes]
        rows = acgh.feature('Row')[valid_probes]
        plt.figure()
        array_image(rows, cols, np.log2(sample_signal), median_center=True)
        plt.figure()
        array_image(rows, cols, np.log2(sample_bg), median_center=True)

        plt.figure()
        array_image(rows, cols, np.log2(reference_signal), median_center=True)
        plt.figure()
        array_image(rows, cols, np.log2(reference_bg), median_center=True)

        plt.figure()
        array_image(rows, cols, M, median_center=True)
        plt.figure()
        array_image(rows, cols, M_norm, median_center=True)

        plt.show()
        exit()

        # Profile plot ----------------------------------------------------------------
        ratio = M_norm

        # Version 1: only ordering of the positions
        plt.figure()
        positions = np.arange(len(ratio))
        separators = np.concatenate(([0], np.cumsum(summary['num_probes'])))
        cgh_profile(positions, ratio[mappings_order], separators=separators)

        # Version 2: ordering and proportional distance between probes
        # Note: the "zero" for each chromosome is the last position of the previous one
        plt.figure()
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

        cgh_profile(positions, ratio, separators=separators)


        # Replicates statistics -------------------------------------------------------
        probe_positions = defaultdict(list)
        for i, probe in enumerate(probes):
            probe_positions[probe].append(i)

        for probe in probe_positions.keys():
            if len(probe_positions[probe]) == 1:
                del probe_positions[probe]

        print
        print 'There are %d replicated probes on the array' % len(probe_positions)

        x = list()
        y = list()
        y_err = list()
        for probe in probe_positions:
            spots = probe_positions[probe]

            x.append(positions[spots][0])
            y.append(ratio[spots].mean())
            y_err.append(ratio[spots].std())

            plt.scatter(positions[spots], ratio[spots], c='k', s=4, edgecolors='none')

        x = np.asarray(x)
        y = np.asarray(y)
        y_err = np.asarray(y_err)
        x_idx = np.argsort(x)

        plt.errorbar(x[x_idx], y[x_idx], y_err[x_idx], fmt=None)

        # Saving csv file -----------------------------------------------------
        # Other info
        #cols = acgh.feature('Col')[valid_probes]
        #rows = acgh.feature('Row')[valid_probes]
        #ratio = np.log2(sample_signal) - np.log2(reference_signal)
        #
        #OUT_PATH = os.path.join(SAMPLES_OUT_DIR, FILE_NAME)
        #with open(OUT_PATH, 'w') as csvout:
        #    header = '\t'.join(["Col", "Row", "Position", "Chromosome",
        #                        "Ref_MedianSignal", "Ref_BGMedianSignal",
        #                        "Sample_MedianSignal", "Sample_BGMedianSignal",
        #                        "LogRatio", "ProbeID"])
        #    csvout.write(header + '\n')
        #    for line in zip(cols, rows, positions, mappings['chromosome'],
        #                    reference_signal, reference_bg,
        #                    sample_signal, sample_bg, ratio, probes):
        #        csvout.write('\t'.join([str(x) for x in line]) + '\n')
        break

plt.show()

