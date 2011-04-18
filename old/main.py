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

    if sample['Sample name'] == 'Sample 35': # 35, 44: spatial bias
    #if sample['platform'] == 'GPL5477':
        # Sample reading ------------------------------------------------------
        # We assume red=ch1=cy5='reference' and green=ch2=cy3='tumor sample'
        swap = False if sample['ch1: source name'] == 'reference' else True
        FILE_NAME = sample['Agilent Feature Extraction file']
        FILE_PATH = os.path.join(SAMPLES_DIR, FILE_NAME)
        acgh = AgilentReader(FILE_PATH)

        print acgh


        # Extact locus and sort by them -----------------------------------------------
        mappings = acgh.toarray(fields=('chromosome', 'start_base', 'end_base'))
        print mappings
        exit()

        mappings_order = np.argsort(mappings) # sort by 'chromosome' and (if equal) by 'start'
                                            # (then, eventually, 'end')

        # Print Sample informations -------------------------------------------------
        print
        print FILE_NAME
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
        plt.figure()
        A = 0.5 * (np.log2(sample_signal) + np.log2(reference_signal))
        M = np.log2(sample_signal) - np.log2(reference_signal)
        M_norm = MA_plot(A, M, lowess=True, label='Raw Signal')

        # Agilent Log Ratio Plot ----------------------------------------------
        plt.figure()
        agilent_M_norm = MA_plot(A, agilent_log_ratio, lowess=True, label='Agilent Signal')

        # Array Image plot ----------------------------------------------------------
        cols = acgh.feature('Col')[valid_probes]
        rows = acgh.feature('Row')[valid_probes]
        #plt.figure()
        #array_image(rows, cols, np.log2(sample_signal), median_center=True)
        #plt.figure()
        #array_image(rows, cols, np.log2(sample_bg), median_center=True)
        #
        #plt.figure()
        #array_image(rows, cols, np.log2(reference_signal), median_center=True)
        #plt.figure()
        #array_image(rows, cols, np.log2(reference_bg), median_center=True)

        #plt.figure()
        #plt.title('Raw Signal')
        #array_image(rows, cols, M, median_center=True)
        plt.figure()
        array_image(rows, cols, M_norm, median_center=True)
        plt.title('Raw Signal (after Lowess)')

        # Agilent Array Image -------------------------------------------------
        plt.figure()
        array_image(rows, cols, agilent_log_ratio, median_center=True)
        plt.title('Agilent Signal')

        # Profile plot ----------------------------------------------------------------
        ratio = M_norm

        # Version 1: only ordering of the positions
        #plt.figure()
        #positions = np.arange(len(ratio))
        #separators = np.concatenate(([0], np.cumsum(summary['num_probes'])))
        #cgh_profile(positions, ratio[mappings_order], separators=separators)

        # Version 2: ordering and proportional distance between probes
        # Note: the "zero" for each chromosome is the last position of the previous one
        #plt.figure()
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

        plt.figure()
        cgh_profile(positions, ratio, separators=separators)
        plt.title('Raw Signal Profile (after Lowess)')

        # Agilent Profile Plot--------------------------------------------------
        plt.figure()
        cgh_profile(positions, agilent_log_ratio, separators=separators)
        plt.title('Agilent Signal Profile')

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

            #plt.scatter(positions[spots], ratio[spots], c='k', s=4, edgecolors='none')

        x = np.asarray(x)
        y = np.asarray(y)
        y_err = np.asarray(y_err)
        x_idx = np.argsort(x)

        #plt.errorbar(x[x_idx], y[x_idx], y_err[x_idx], fmt=None)

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

        from manor import manor_normalization
        arrayCGH_norm = manor_normalization(acgh, valid_probes, probes, reference_signal,
                            reference_bg, sample_signal, sample_bg, M,
                            mappings, positions)

        new_valid_probes = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Flag')) == 'OK'
        ratio_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('LogRatioNorm'))[new_valid_probes]
        rows = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Row'))[new_valid_probes]
        cols = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Col'))[new_valid_probes]
        sample_signal_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Sample_MedianSignal'))[new_valid_probes]
        reference_signal_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Ref_MedianSignal'))[new_valid_probes]

        #plt.figure()
        #plt.title('Spatial Normalized Signal')
        #array_image(rows, cols, ratio_norm, median_center=True)

        plt.figure()
        A = 0.5 * (np.log2(sample_signal_norm) + np.log2(reference_signal_norm))
        M = np.log2(sample_signal_norm) - np.log2(reference_signal_norm)
        M_s_norm = MA_plot(A, M, lowess=True, label='Spatial Normalized Signal')

        plt.figure()
        array_image(rows, cols, M_s_norm, median_center=True)
        plt.title('Spatial filtered probes (after Lowess)')

        positions_norm = np.empty(new_valid_probes.shape)
        positions_norm[valid_probes] = positions
        positions_norm = positions_norm[new_valid_probes]

        plt.figure()
        cgh_profile(positions_norm, M_s_norm, separators=separators)
        plt.title('Spatial filtered Signal Profile (after Lowess)')

        # ---------------------------------------------------------------------

        from rpy2 import robjects
        from rpy2.robjects.packages import importr
        importr('CGHnormaliter')

        rvalid_probes = robjects.BoolVector(new_valid_probes)
        probes_id = arrayCGH_norm.rx2('arrayValues').rx2('ProbeID').rx(rvalid_probes)
        chromosomes = arrayCGH_norm.rx2('arrayValues').rx2('Chromosome').rx(rvalid_probes)

        starts = np.empty(new_valid_probes.shape)
        starts[valid_probes] = mappings['start']
        starts = starts[new_valid_probes]

        ends = np.empty(new_valid_probes.shape)
        ends[valid_probes] = mappings['end']
        ends = ends[new_valid_probes]

        probes_norm = np.asarray(np.empty(new_valid_probes.shape), dtype=probes.dtype)
        probes_norm[valid_probes] = probes
        probes_norm = probes_norm[new_valid_probes]

        input_names = ['CloneID', 'Chromosome', 'Start', 'End', 'Case1.test',
                       'Case1.ref']
        unique_probes, unique_index = np.unique(probes_norm, return_index=True)
        runique_index = robjects.IntVector(unique_index+1) #R starts from 1
        input_values = [probes_id.rx(runique_index),
                        robjects.IntVector(np.asarray(chromosomes)).rx(runique_index),
                        robjects.IntVector(starts).rx(runique_index),
                        robjects.IntVector(ends).rx(runique_index),
                        robjects.FloatVector(sample_signal_norm).rx(runique_index),
                        robjects.FloatVector(reference_signal_norm).rx(runique_index)]

        import rpy2.rlike.container as rlc
        d = rlc.OrdDict(zip(input_names, input_values))
        normaliter_data = robjects.DataFrame(d)

        CGHnormaliter = robjects.r['CGHnormaliter']

        #rwrite_table = robjects.r['write.table']
        #rwrite_table(normaliter_data, "myFile.csv", col_names = True,
                     #sep = "\t", row_names=False, quote=False)

        result = CGHnormaliter(normaliter_data)
        row_names = list(robjects.r['copynumber'](result).rownames)
        normalized = np.asarray(robjects.r['copynumber'](result)).squeeze()
        segmented = np.asarray(robjects.r['segmented'](result)).squeeze()
        called = np.asarray(robjects.r['calls'](result)).squeeze()

        #ordered_positions = dict(zip(probes, positions))
        #new_positions = list()
        #for id in row_names:
        #    new_positions.append(ordered_positions[id])

        plt.figure()
        cgh_profile(positions_norm[unique_index], normalized, separators=separators)
        plt.title('Signal Profile (after CGHNormaliter)')
        #plt.scatter(positions_norm[unique_index], called)


        plt.figure()
        order = np.argsort(positions_norm[unique_index])
        cgh_profile(positions_norm[unique_index][order],
                    (normalized - M_s_norm[unique_index])[order],
                    separators=separators)
        plt.title('Difference Signal Profile (CGHNormaliter - Lowess)')


        # Distributions
#plt.hist(M_norm, 100, normed=1, alpha=0.7, label='Raw Lowess')
#plt.hist(M_s_norm, 100, normed=1, alpha=0.7, label='Spatial Filtering + Lowess')
plt.hist(normalized, 100, normed=1, alpha=0.7, label='Spatial Filtering + CGHNormaliter')
plt.hist(agilent_log_ratio, 100, normed=1, alpha=0.7, label='Agilent')
plt.legend()
#plt.axvline(0.0, color='k', lw=2)


plt.show()
