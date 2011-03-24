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
CLINICAL_INFO_PATH = os.path.join(ROOT_DIR, 'info_Row_aCGH.csv')


# As example we pick the first GPL5477 sample reading the clinical data -------
clinical_info = DictReader(open(CLINICAL_INFO_PATH, 'rb'), dialect='excel')
for sample in clinical_info:
    if sample['platform'] == 'GPL5477': break
    #if sample['Sample name'] == 'Sample 3': break

# Sample reading --------------------------------------------------------------
# We assume red=ch1=cy5='reference' and green=ch2=cy3='tumor sample'
swap = False if sample['ch1: source name'] == 'reference' else True
FILE_NAME = sample['Agilent Feature Extraction file']
FILE_PATH = os.path.join(SAMPLES_DIR, FILE_NAME)
acgh = AgilentReader(FILE_PATH)

# Remove controls and unmapped chr and extract data ---------------------------
valid_probes = probes_filter(acgh)
locations = acgh.feature('SystematicName')[valid_probes]
probes = acgh.feature('ProbeName')[valid_probes]

# Default: sample on the green=ch2=cy3 channel
ch1_signal = acgh.feature('rMedianSignal')[valid_probes]
ch2_signal = acgh.feature('gMedianSignal')[valid_probes]
if sample['ch1: source name'] == 'reference':
    reference_signal, sample_signal = ch1_signal, ch2_signal
else:
    reference_signal, sample_signal = ch2_signal, ch1_signal

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
plt.figure()
A = 0.5 * (np.log2(sample_signal) + np.log2(reference_signal))
M = np.log2(sample_signal) - np.log2(reference_signal)
M_norm = MA_plot(A, M, lowess=True, label='Raw Signal')

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

cgh_profile(positions, ratio, separators=separators)

plt.show()


exit()

# Replicates statistics -------------------------------------------------------
probe_positions = defaultdict(list)
for i, probe in enumerate(probes):
    probe_positions[probe].append(i)

for probe in probe_positions.keys():
    if len(probe_positions[probe]) == 1:
        del probe_positions[probe]

print
print 'There are %d replicated probes on the array' % len(probe_positions)

# Neuvial et Al. aCGH data quality assessment
# Sigma: experimental noise (median of the std. devs log ratio)
stds = list()
for probe in probe_positions:
    pos = probe_positions[probe]
    stds.append(ratio[pos].std(ddof=1))

sigma = np.median(stds)
print 'aCGH sigma (agilent log ratios): %f' % sigma

# Spatial plot ---------------------------------------------------------------
rows = acgh.feature('Row')[valid_probes]
cols = acgh.feature('Col')[valid_probes]
#fig = array_image(rows, cols, ratio, num=1)


## Proviamo con la griglia GEO
#agilent_feat_id = acgh.feature('FeatureNum')[valid_probes]
#gpl_feat_id = gpl.field('ID')
#gpl_cols = gpl.field('COL')
#gpl_rows = gpl.field('ROW')
#
#num_rows = gpl_rows.max()
#num_cols = gpl_cols.max()
#
#num_ratio_mapping = dict()
#for fn, r in zip(agilent_feat_id, ratio):
#    num_ratio_mapping[fn] = r
#
#gpl_array_image = np.ones((num_rows, num_cols)) * np.inf
#for r, c, id in zip(gpl_rows, gpl_cols, gpl_feat_id):
#    gpl_array_image[r-1, c-1] = num_ratio_mapping.get(id, np.inf)


# M-A plots -------------------------------------------------------------------
#plt.figure()
#plt.title('Agilent normalization')
#A = 0.5 * (np.log2(red) + np.log2(green))
#M = np.log2(red) - np.log2(green) # == ratio
#plt.scatter(A, M, s=4)
#plt.axhline(np.median(M))
#
#plt.figure()
#plt.title('Median Signal')
#green = acgh.feature('gMedianSignal')[valid_probes]
#red = acgh.feature('rMedianSignal')[valid_probes]
#A = 0.5 * (np.log2(red) + np.log2(green))
#M = np.log2(red) - np.log2(green)
#plt.scatter(A, M, s=4)
#plt.axhline(np.median(M))
#
#plt.figure()
#plt.title('BG Subtracted Signal')
#green = acgh.feature('gBGSubSignal')[valid_probes]
#red = acgh.feature('rBGSubSignal')[valid_probes]
#A = 0.5 * (np.log2(red) + np.log2(green))
#M = np.log2(red) - np.log2(green)
#plt.scatter(A, M, s=4)
#plt.axhline(np.median(M))
#
#plt.show()

# TODO
# Plot:
    #- M-A
    #- Vedi GenomeGraphs in R (curtis)
    #- spatial Plots (Neuvial, Khojasteh) ()
    #- Log2 colorato (//)

# Metriche
    #- Sigma, smt, dyn (Neuvial)






#for chr, start, end in loc_tuples:
    #print chr, end-start


## Plots mappings
#plt.figure()
#
#xmin = np.inf
#xmax = -np.inf
#for ch, start, end in loc_tuples:
#    if ch == 'chr1':
#        plt.plot([int(start), int(end)], [1.0, 1.0], '.', lw=5)
#        plt.plot([int(start), int(end)], [1.0, 1.0], '-', lw=2)
#        if int(start) < xmin: xmin = int(start)
#        if int(end) > xmax: xmax = int(end)
#
##plt.plot([xmin, xmin+10], [1.0, 1.0], 'k-', lw=1)
#
##print xmin, xmax
#plt.axis([xmin, xmax, 0.5, 1.5])
#
##ax.set_xticks([])
##ax.set_yticks([])
#plt.show()


#
## Plot Boxplots of duplicated probes
#duplicated = list()
#for c in counter:
#    if len(counter[c]) > 1:
#        duplicated.append(counter[c])
#
#plt.boxplot(duplicated)
