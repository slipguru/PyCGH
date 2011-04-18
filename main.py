import os
from collections import defaultdict
from csv import DictReader

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab as ml

from cghutils.readers import AgilentReader, GPLReader
from cghutils.filters import probes_filter, split_mappings
from cghutils.plots import array_image, MA_plot, cgh_profile

ROOT_DIR = '/home/sabba/Phd/Tonini_IST'
SAMPLES_DIR = os.path.join(ROOT_DIR, 'Group1')
SAMPLES_OUT_DIR = os.path.join(ROOT_DIR, 'Group1_extracted')
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')


out_replicated = list()
out_names = list()

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
        acgh_r = AgilentReader(FILE_PATH)
        acgh_full = acgh_r.toarray(order=('row', 'col')) # def
        acgh = acgh_full[acgh_full['valid']]

        #=======================================
        #from mlabwrap import mlab
        #mlab.addpath('chen/')
        #mlab.addpath('chen/stprtool/')
        #mlab.stprpath('chen/stprtool/')
        #data = mlab.ReadAgilentResult(FILE_PATH)
        #aCGHdata, model = mlab.aCGHNormalization3(data, 1, (0 if swap else 1), 1, nout=2)
        #plt.scatter(np.arange(len(aCGHdata.nRatioAdj)),
        #            aCGHdata.nRatioAdj.flatten(),
        #            c=aCGHdata.nRatioAdj.flatten(),
        #            cmap=plt.get_cmap('jet'),
        #            vmin=-1, vmax=1, s=8, edgecolors='none')
        #plt.show()
        #exit()
        #=======================================


        if swap:
            test_signal, reference_signal = acgh['ref_signal'], acgh['test_signal']
        else:
            test_signal, reference_signal = acgh['test_signal'], acgh['ref_signal']
        M = np.log2(test_signal) - np.log2(reference_signal)
        plt.figure()
        array_image(acgh['row'], acgh['col'], M, median_center=True)

        #plt.figure()
        #A = reference_signal
        #M = test_signal
        #lws_ratio = MA_plot(A, M, lowess=True, label='Raw Signal')

        plt.show()
        exit()

        #======================================================================
        #if swap:
        #    test_signal, reference_signal = acgh['ref_signal'], acgh['test_signal']
        #else:
        #    test_signal, reference_signal = acgh['test_signal'], acgh['ref_signal']
        #
        #A = 0.5 * (np.log2(test_signal) + np.log2(reference_signal))
        #M = np.log2(test_signal) - np.log2(reference_signal)
        #lws_ratio = MA_plot(A, M, lowess=True, label='Raw Signal')
        #
        ## Replicates statistics -------------------------------------------------------
        #probes = acgh['id']
        #probe_positions = defaultdict(list)
        #for i, probe in enumerate(probes):
        #    probe_positions[probe].append(i)
        #
        #for probe in probe_positions.keys():
        #    if len(probe_positions[probe]) == 1:
        #        del probe_positions[probe]
        #
        #print
        #print 'There are %d replicated probes on the array' % len(probe_positions)
        #
        #values = list()
        #for probe in probe_positions:
        #    spots = probe_positions[probe]
        #    values.append(lws_ratio[spots].std()) #deviazione std dei replicati
        #
        #out_replicated.append(np.array(values))
        #out_names.append(int(sample['Sample name'].split()[1]))
        #continue
        #======================================================================

        # Print Sample informations -------------------------------------------
        print
        print FILE_NAME
        print '\n'.join(('Sample name: %(Sample name)s',
                         'Platform: %(platform)s')) % sample

        # Mappings fields are: 'chrmosome', 'start_base', 'end_base'
        summary = ml.rec_groupby(acgh_full,
                                   groupby=('chromosome',),
                                   stats=(('chromosome', len, 'num_probes'),
                                          ('start_base', np.min, 'from'),
                                          ('end_base', np.max, 'to'),
                                          ))
        rjusts = (6, 10, 10, 10)
        print 'Chr # | # probes |   from   |    to    '
        print '------+----------+----------+----------'
        for record in summary:
            formatted = "|".join(str(val).rjust(just) for val, just in zip(record, rjusts))
            print formatted
        summary = summary[1:]

        # Chromosomes position on the glass -----------------------------------
        array_image(acgh['row'], acgh['col'], acgh['chromosome'],
                    median_center=False, vmin=1, vmax=24)
        plt.title('Chromosome position on the glass')
        #plt.show()
        #exit()

        # Signals extraction --------------------------------------------------
        if swap:
            test_signal, reference_signal = acgh['ref_signal'], acgh['test_signal']
        else:
            test_signal, reference_signal = acgh['test_signal'], acgh['ref_signal']


        #plt.figure()
        #M = np.log2(test_signal) - np.log2(reference_signal)
        #plt.hist(normalized, 100, normed=1, alpha=0.7, label='Spatial Filtering + CGHNormaliter')
        #plt.hist(M, 100, normed=False, alpha=0.7, label='Raw')
        #plt.legend()
        #plt.show()
        #continue

        #exit()

        # M-A plot ------------------------------------------------------------
        plt.figure()
        A = 0.5 * (np.log2(test_signal) + np.log2(reference_signal))
        M = np.log2(test_signal) - np.log2(reference_signal)
        lws_ratio = MA_plot(A, M, lowess=True, label='Raw Signal')

        plt.figure()
        plt.plot(test_signal, reference_signal, '.')

        # Agilent Log Ratio Plot ----------------------------------------------
        #plt.figure()
        #agilent_M_norm = MA_plot(A, agilent_log_ratio, lowess=True, label='Agilent Signal')

        # Array Image plot ----------------------------------------------------------
        plt.figure()
        array_image(acgh['row'], acgh['col'], lws_ratio, median_center=True)
        plt.title('Raw Signal (after Lowess)')

        # Agilent Array Image -------------------------------------------------
        #plt.figure()
        #array_image(rows, cols, agilent_log_ratio, median_center=True)
        #plt.title('Agilent Signal')

        # Profile plot ----------------------------------------------------------------
        # Ordering and proportional distance between probes
        # Note: the "zero" for each chromosome is the last position of the previous one
        plt.figure()
        positions = np.empty_like(lws_ratio)
        separators = np.concatenate(([0], np.cumsum(summary['to'])))

        # Shifting...
        mappings = acgh[['chromosome', 'start_base', 'end_base']]
        starting_locations = dict(zip(summary['chromosome'], separators[:-1]))
        for i, (chr, start, end) in enumerate(mappings):
            shifted_start = start + starting_locations[chr]
            shifted_end = end + starting_locations[chr]
            pos = int((shifted_start + shifted_end)/2) # median point of the probe
            positions[i] = pos

        # PROBLEMA: le probe replicate hanno stessa position...
        # lo scatter plot li stampa tutti (dorvrebbero essere visibili)
        plt.figure()
        cgh_profile(positions, lws_ratio, separators=separators)
        plt.title('Raw Signal Profile (after Lowess)')

        srtindexes = np.argsort(mappings)
        plt.figure()
        cgh_profile(positions[srtindexes][:len(aCGHdata.nRatioAdj.flatten())],
                    aCGHdata.nRatioAdj.flatten(),
                    separators=separators)#, vmin=-.6, vmax=.6)
        plt.show()
        exit()

        # Agilent Profile Plot--------------------------------------------------
        #plt.figure()
        #cgh_profile(positions, agilent_log_ratio, separators=separators)
        #plt.title('Agilent Signal Profile')

        # Replicates statistics -------------------------------------------------------
        probes = acgh['id']
        probe_positions = defaultdict(list)
        for i, probe in enumerate(probes):
            probe_positions[probe].append(i)

        for probe in probe_positions.keys():
            if len(probe_positions[probe]) == 1:
                del probe_positions[probe]

        print
        print 'There are %d replicated probes on the array' % len(probe_positions)

        replicate_indexed = np.ones_like(lws_ratio) * np.nan

        x = list()
        y = list()
        y_err = list()
        for i, probe in enumerate(probe_positions):
            spots = probe_positions[probe]

            x.append(positions[spots][0])
            y.append(lws_ratio[spots].mean())
            #y_err.append(lws_ratio[spots].std())
            if i == 0:
                replicate_indexed[spots] = (i+1)

            #plt.scatter(positions[spots], lws_ratio[spots], c='k', s=8, edgecolors='none')

        x = np.asarray(x)
        y = np.asarray(y)
        #y_err = np.asarray(y_err)
        x_idx = np.argsort(x)
        #plt.errorbar(x[x_idx], y[x_idx], y_err[x_idx], fmt=None)
        #plt.scatter(x[x_idx], y[x_idx], c='r', s=8)
        #print replicate_indexed

        plt.figure()
        array_image(acgh['row'], acgh['col'], replicate_indexed,
                    median_center=False, vmin=1, vmax=(i+1))
        plt.title('Replicate positions')
        #exit()


        #from scikits.learn import mixture
        #clf = mixture.GMM(n_states=2, cvtype='full')
        #X = np.c_[M]#, A] # dati sono [A, M]
        #clf.fit(X)
        #Y = clf.predict(X)
        #means = clf.means.flatten()
        #covars = clf.covars.flatten()
        #points = np.linspace(M.min(), M.max(), 1000)
        #plt.plot(points, plt.normpdf(points, means[0], covars[0]))
        #plt.plot(points, plt.normpdf(points, means[1], covars[1]))
        #
        #plt.figure()
        #plt.scatter(A[Y==0], M[Y==0], s=8, c='r', edgecolors='none', alpha=.7)
        #plt.scatter(A[Y==1], M[Y==1], s=8, c='b', edgecolors='none', alpha=.7)
        #normals = 0 if (Y==0).sum() > (Y==1).sum() else 1
        #print normals
        #
        #plt.figure()
        #lws_ratio_ = MA_plot(A[Y==normals], M[Y==normals], lowess=True, label='Raw Signal')
        #lws_ratio2 = lws_ratio.copy()
        #lws_ratio2[Y==normals] = lws_ratio_
        #
        #plt.figure()
        #cgh_profile(positions, lws_ratio2, separators=separators)

        #plt.show()
        #exit()

        # MANOR ---------------------------------------------------------------
        from manor import manor_normalization

        arrayCGH_norm = manor_normalization(acgh_full, M, positions)

        new_valid_probes = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Flag')) == 'OK'
        #new_valid_probes = np.logical_and(new_valid_probes, acgh_full['valid'])
        ratio_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('LogRatioNorm'))[new_valid_probes]
        rows = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Row'))[new_valid_probes]
        cols = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Col'))[new_valid_probes]
        test_signal_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Sample_MedianSignal'))[new_valid_probes]
        reference_signal_norm = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2('Ref_MedianSignal'))[new_valid_probes]
        trend = np.asarray(arrayCGH_norm.rx2('arrayValues').rx2("Trend"))#[new_valid_probes]

        plt.figure()
        plt.title('Spatial Normalized Signal')
        array_image(rows, cols, ratio_norm, median_center=True)

        plt.figure()
        array_image(rows, cols, trend, median_center=True, vmin=-0.3, vmax=.3)

        #plt.figure()
        A = 0.5 * (np.log2(test_signal_norm) + np.log2(reference_signal_norm))
        M = np.log2(test_signal_norm) - np.log2(reference_signal_norm)
        M_s_norm = MA_plot(A, M, lowess=True, label='Spatial Normalized Signal')

        plt.figure()
        array_image(rows, cols, ratio_norm, median_center=True)
        plt.title('Spatial filtered probes (after Lowess)')

        positions_norm = np.empty(new_valid_probes.shape)
        positions_norm[acgh_full['valid']] = positions
        positions_norm = positions_norm[new_valid_probes]

        plt.figure()
        cgh_profile(positions_norm, M, separators=separators)
        plt.title('Spatial filtered Signal Profile (after Lowess)')

        # ---------------------------------------------------------------------

        from rpy2 import robjects
        from rpy2.robjects.packages import importr
        importr('CGHnormaliter')

        rvalid_probes = robjects.BoolVector(new_valid_probes)
        probes_id = arrayCGH_norm.rx2('arrayValues').rx2('ProbeID').rx(rvalid_probes)
        chromosomes = arrayCGH_norm.rx2('arrayValues').rx2('Chromosome').rx(rvalid_probes)

        starts = np.empty(new_valid_probes.shape)
        starts[acgh_full['valid']] = mappings['start_base']
        starts = starts[new_valid_probes]

        ends = np.empty(new_valid_probes.shape)
        ends[acgh_full['valid']] = mappings['end_base']
        ends = ends[new_valid_probes]

        probes_norm = np.asarray(np.empty(new_valid_probes.shape), dtype=probes.dtype)
        probes_norm[acgh_full['valid']] = probes
        probes_norm = probes_norm[new_valid_probes]

        input_names = ['CloneID', 'Chromosome', 'Start', 'End', 'Case1.test',
                       'Case1.ref']
        unique_probes, unique_index = np.unique(probes_norm, return_index=True)
        runique_index = robjects.IntVector(unique_index+1) #R starts from 1
        input_values = [probes_id.rx(runique_index),
                        robjects.IntVector(np.asarray(chromosomes)).rx(runique_index),
                        robjects.IntVector(starts).rx(runique_index),
                        robjects.IntVector(ends).rx(runique_index),
                        robjects.FloatVector(test_signal_norm).rx(runique_index),
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


plt.figure()
plt.boxplot(out_replicated, notch=False, positions=range(len(out_names)))
plt.xticks(range(len(out_names)), out_names)
plt.show()
