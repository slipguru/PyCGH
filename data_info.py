import os
from csv import DictReader

import numpy as np
from matplotlib import pyplot as plt

from cghreader import CGHReader

DATA_PATH = './aCGH/'

def _main():
    info_csv = DictReader(open('info_Row_aCGH.csv', 'rb'), dialect='excel')

    for sample in info_csv:
        aCGH_path = os.path.join(DATA_PATH,
                                 sample['Agilent Feature Extraction file'])

        # Sample corrotto
        if sample['Sample name'] == 'Sample 21':
            continue
        #if sample['platform'] != 'GPL2873':
            #continue

        #if sample['Sample name'] != 'Sample 3':
            #continue

        acgh = CGHReader(aCGH_path)
        fe_version = acgh.FEPARAMS['FeatureExtractor_Version']
        time = acgh.FEPARAMS['FeatureExtractor_ExtractionTime']
        num_tot_feat = len(acgh.FEATURES['FeatureNum'])
        controls = (acgh.FEATURES['ControlType'] == 1)
        num_feat = num_tot_feat - len(acgh.FEATURES['FeatureNum'][controls])
        num_unmapped = (acgh.FEATURES['SystematicName'] == 'unmapped').sum()


        filename = sample['Agilent Feature Extraction file']
        print '\t'.join((sample['Sample name'],
                         sample['platform'],
                         fe_version,
                         time,
                         str(num_tot_feat),
                         str(num_feat),
                         filename
                         #str(num_unmapped)
                        ))

    #def print_log_ratio(acgh):
        # Print LogRatio test
        syst_names = acgh.FEATURES['SystematicName']
        controls = acgh.FEATURES['ControlType']
        mapped = np.array([i for i, v in enumerate(syst_names)
                                         if v.startswith('ch') and
                                         not 'random' in v and
                                         not 'hla_hap1' in v and
                                         controls[i] == 0], dtype=int)

        syst_names = syst_names[mapped]
        log_ratio = acgh.FEATURES['LogRatio'][mapped]
        #plt.plot(log_ratio, 'r.')
        #plt.figure()

        def extract_ch_key(pair):
            value = pair[1]
            ch, start = value.split(':')
            ch = ch.replace('chr', '')
            try:
                ch = int(ch)
            except ValueError:
                pass

            start = int(start.split('-')[0])

            return (ch, start)

        zipped_solution = sorted(zip(log_ratio, syst_names), key=extract_ch_key)
        sorted_log_ratio, sorted_chr = zip(*zipped_solution)

        # Aree cromosomi
        limits = dict()
        for i, chr in enumerate(sorted_chr):
            ch = chr.split(':')[0].replace('chr', '')
            try:
                ch = int(ch)
            except ValueError:
                pass

            limits[ch] = i


        # Plot tutti cromosomi -----------------------------
        plt.figure(1)
        plt.plot(sorted_log_ratio, 'b.')
        plt.axhline(linewidth=2, color='r', alpha=0.5)
        plt.axhline(y=1, linewidth=1, color='r', alpha=0.5)
        plt.axhline(y=-1, linewidth=1, color='r', alpha=0.5)
        plt.ylim(-1.5, 1.5)

        y_annot = log_ratio.min() - 0.3
        prev = 0
        ticks = []
        for l in sorted(limits.keys()):
            plt.axvline(limits[l], color='k')

            ticks.append(prev + ((limits[l]-prev)/2))
            prev = limits[l]

        plt.xticks(ticks,
                   ['Chr %s' % k for k in sorted(limits.keys())],
                   rotation=90)

        # Plots per singolo cromosoma -------------------------
        plt.figure(2)
        prev = 0
        for i, l in enumerate(sorted(limits.keys())):
            plt.subplot(6, 4, i+1)
            plt.title(l)

            plt.plot(sorted_log_ratio[prev:limits[l]], 'b.')
            plt.axhline(linewidth=1, color='r', alpha=0.5)
            plt.axhline(y=1, linewidth=0.5, color='r', alpha=0.5)
            plt.axhline(y=-1, linewidth=0.5, color='r', alpha=0.5)
            plt.ylim(-1.5, 1.5)

            prev = limits[l]

        plt.show()


if __name__ == '__main__':
    main()

#./aCGH/18086_G1_Cologne_251469812734_1_1.txt
