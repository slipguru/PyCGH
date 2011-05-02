#-*- coding: utf-8 -*-

import os

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from matplotlib import pylab as plt
from cghutils.plots import array_image, MA_plot, cgh_profile
from cghutils import readers

# Importazione dati
importr('MANOR')
robjects.r('data(spatial)')

def manor_normalization(acgh, ratio, positions):
    # ArrayValues -----------------------------------------------------------------
    spot_names = ["Col", "Row", "ProbeID",
                  "Ref_MedianSignal", "Ref_BGMedianSignal",
                  "Sample_MedianSignal", "Sample_BGMedianSignal",
                  "LogRatio", "Chromosome"]

    # Conversione, per gestire anche i valori missing NA richiesti da MANOR
    def _convert(indexes, values, default, function):
        conv = np.asarray(np.empty_like(indexes), dtype=object)
        conv[:] = default
        conv[indexes] = values
        return function(conv)

    def convert_string(indexes, values):
        return _convert(indexes, values, robjects.NA_character, robjects.StrVector)

    def convert_float(indexes, values):
        return _convert(indexes, values, robjects.NA_real, robjects.FloatVector)

    def convert_int(indexes, values):
        return _convert(indexes, values, robjects.NA_integer, robjects.IntVector)

    cols = acgh['col']#.feature('Col')
    rows = acgh['row']#.feature('Row')
    cols_max = cols.max()
    rows_max = rows.max()

    spot_values = [robjects.IntVector(cols),
                   robjects.IntVector(rows),
                   convert_string(acgh['valid'], acgh['id']),
                   convert_float(acgh['valid'], acgh['ref_signal']),
                   convert_float(acgh['valid'], acgh['ref_bg']),
                   convert_float(acgh['valid'], acgh['test_signal']),
                   convert_float(acgh['valid'], acgh['test_bg']),
                   convert_float(acgh['valid'], ratio),
                   convert_int(acgh['valid'], acgh['chromosome'])]

    arrayValues = robjects.DataFrame(dict(zip(spot_names, spot_values)))

    for k in spot_names:
        print k, len(arrayValues.rx2(k))

    # CloneValues -----------------------------------------------------------------
    probes = acgh['id'][acgh['valid']]
    chromosomes = acgh['chromosome'][acgh['valid']]
    clone_names = ["ProbeID", "Position", "Chromosome"]
    unique_probes, unique_index = np.unique(probes, return_index=True)
    clone_values = [robjects.FactorVector(unique_probes),
                    robjects.IntVector(positions[unique_index]),
                    robjects.IntVector(chromosomes[unique_index])]
    cloneValues = robjects.DataFrame(dict(zip(clone_names, clone_values)))

    for k in clone_names:
        print k, len(cloneValues.rx2(k))

    # ArrayCgh --------------------------------------------------------------------
    rlist = robjects.r('list')
    arrayCGH = rlist(**{'arrayValues': arrayValues,
                        'cloneValues': cloneValues,
                        'arrayDesign': robjects.IntVector([1, 1, cols_max, rows_max]),
                        'id.rep': "ProbeID"
                        })
    robjects.globalenv['acgh'] = arrayCGH # put the arrayCGH in the R environment
    robjects.r('class(acgh) <- "arrayCGH"')
    arrayCGH = robjects.r['acgh']

    # Normalization ---------------------------------------------------------------
    robjects.r("source('cghutils/manor_normalization.R')")

    manor_norm = robjects.r['manor_norm']
    arrayCGH_norm = manor_norm(arrayCGH)

    return arrayCGH_norm
