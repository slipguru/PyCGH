#-*- coding: utf-8 -*-

import os

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc

from matplotlib import pylab as plt
from cghutils.plots import array_image, MA_plot, cgh_profile

def manor_test(ratio, ref_fore, ref_bg, sample_fore, sample_bg, pos_order,
               chromosomes):
    importr('MANOR')
    print robjects.r('sessionInfo()')

    # Importazione dati
    robjects.r('data(flags)')
    robjects.r('data(spatial)')

    # Valori richiesti nella lista:
    # arrayValues
    #       LogRatio, RefFore, RefBack, Dapifore, DapiBack, PosOrder
    # cloneValues
    #       PosOrder, Chromosome
    # id.rep

    od1 = rlc.OrdDict([('LogRatio', ratio),
                      ('RefFore', ref_fore),
                      ('RefBack', ref_bg),
                      ('DapiFore', sample_fore),
                      ('DapiBack', sample_bg),
                      ('PosOrder', pos_order)])
    arrayValues = robjects.DataFrame(od1)
    od2 = rlc.OrdDict([('PosOrder', pos_order),
                      ('Chromosome', chromosomes),
                      ('RefBack', ref_bg)])
    arrayValues = robjects.DataFrame(od1)



    # Questi li abbiamo o possiamo calcolarli (dubbi su ScaledLogRatio ??)
    spot_names = np.asarray(["LogRatio", "RefFore", "RefBack", "DapiFore",
                             "DapiBack", "SpotFlag", "ScaledLogRatio"])

    # PosOrder indica l'ordinamento delle probe e PosOrder uguali indicano
    # probes replicate (non serve posizionamento fisico, ma può essere utilizzato
    # come pos order anche se non contigui)
    # Quindi PosOrder è utilizzato come identificativo del clone
    clone_names = np.asarray(["PosOrder", "Chromosome"])

    # Importazione dati esempio
    acgh = robjects.r['edge']

    # La funzione ritorna una lista di 4 elementi
    # edge$arrayValues  edge$arrayDesign  edge$cloneValues  edge$id.rep
    print acgh.names
    print acgh[0].names
    plt.figure()
    array_image(np.asarray(acgh[0][1]), np.asarray(acgh[0][0]), np.asarray(acgh[0][2]))
    #plt.show()
    #exit()

    # Configurazione flag normalizzaione spaziale
    robjects.r('local.spatial.flag$args <- alist('
               'var="ScaledLogRatio", '     # Variabile da normalizzare
               'by.var=NULL, '              # Centratura rispetto a questa di var
               'nk=5, '                     # Nem: numbero of classes
               'prop=0.25, '                # Massima proporzione di array che
                                            # può essere affetta da bias spaziale
               'thr=0.15, '                 # Thr per spatial bias
               'beta=1, '                   # Nem: Scale parameter
               'family="gaussian")')        # Loess (gaussian -> least-squares)
    robjects.r('flag.list <- list('
               'spatial=local.spatial.flag, '   # Normalizzazione spaziale (S) "ESCLUSIONE"
               #'spot=spot.corr.flag, '          # Negative spot correlation (C)
               'ref.snr = ref.snr.flag, '       # Low S/N (Ref) con thr=1.25 [RefFore/RefBack] (B)
               'dapi.snr = dapi.snr.flag, '     # Low S/N (Dapi) con thr=1.25 [DapiFore/DapiBack] (D)
               'rep=rep.flag, '                 # Poor replicate consistency thr=0.1 (E)
               'unique=unique.flag)')           # Unique Spot (U)
    flag_list = robjects.r['flag.list']

    # Funzione di normalizzazione
    norm = robjects.r['norm']
    median = robjects.r['median']
    acgh_norm = norm(acgh, **{'flag.list':flag_list, 'FUN': median, 'na.rm':True})
    flags = np.asarray(acgh_norm[0][11]) #flags
    not_ok = np.logical_not(flags == 'OK')

    print acgh_norm.names
    print acgh_norm[0].names
    plt.figure()
    array_image(np.asarray(acgh_norm[0][1]), np.asarray(acgh_norm[0][0]), np.asarray(acgh_norm[0][10]))
    print np.asarray(acgh_norm[0][10])[not_ok]
    print np.asarray(acgh[0][9])[not_ok]

    before_pos = np.asarray(acgh[0][9]) # arrayValues$PosOrder
    before_ratios = np.asarray(acgh[0][10]) # arrayValues$LogRatio

    after_pos = np.asarray(acgh_norm[0][9])
    after_ratios = np.asarray(acgh_norm[0][10]) # LogRatioNorm

    plt.figure()
    cgh_profile(before_pos, before_ratios, np.asarray([before_pos.min(), before_pos.max()]))

    plt.figure()
    cgh_profile(after_pos, after_ratios, np.asarray([after_pos.min(), after_pos.max()]))

    plt.show()


    #print before.names

    #edge.norm <- norm(edge, flag.list=flag.list, FUN = median, na.rm = TRUE)
