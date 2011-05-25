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
