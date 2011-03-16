import os
from collections import defaultdict

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab

from cghutils.readers import AgilentReader, GPLReader
from cghutils.filters import probes_filter, split_locus_mapping
from cghutils.plots import array_image, MA_plot, cgh_profile

ABS_DIR = '/home/sabba/Phd/Tonini_IST/Group1/'
FILE_NAME = '1115_G1_251495014936_1_2.txt'          # GPL5477
#FILE_NAME = '1431_G1_251495014937_1_1.txt'
#FILE_NAME = '1912_G1_251495014704_1_2.txt'
#FILE_NAME = '2024_G1_251495014939_1_2.txt'
#FILE_NAME = '2216_G1_251495014704_1_3.txt'
#FILE_NAME = '2717_G1_251495014464_1_2.txt'
#FILE_NAME = '17869_G1_Cologne_251469812733_1_2.txt' # GPL 4093
ABS_PATH = os.path.join(ABS_DIR, FILE_NAME)

acgh = AgilentReader(ABS_PATH)

# GEO GPL reading
#gpl = GPLReader('/home/sabba/Phd/Tonini_IST/Platforms/GPL5477.txt')
#gpl = GPLReader('/home/sabba/Phd/Tonini_IST/Platforms/GPL4093.txt')

# Remove controls and unmapped chr and extract data ---------------------------
valid_probes = probes_filter(acgh)
locations = acgh.feature('SystematicName')[valid_probes]
probes = acgh.feature('ProbeName')[valid_probes]

#green = acgh.feature('gMedianSignal')[valid_probes]
#red = acgh.feature('rMedianSignal')[valid_probes]
#green = acgh.feature('gBGSubSignal')[valid_probes]
#red = acgh.feature('rBGSubSignal')[valid_probes]
red = acgh.feature('rProcessedSignal')[valid_probes]
green = acgh.feature('gProcessedSignal')[valid_probes]
ratio = np.log2(red/green)

# Extact locus and sort by them -----------------------------------------------
locus = split_locus_mapping(locations) # locus is a record array
locus_order = np.argsort(locus) # sort by 'chromosome' and (if equal) by 'start'
                                # (then, eventually, 'end')

# Print a chromosome summary of the start end interval ------------------------
# Locus fields are: 'chrmosome', 'start', 'end'
summary = mlab.rec_groupby(locus,
                           groupby=('chromosome',),
                           stats=(('start', np.min, 'from'),
                                  ('end', np.max, 'to'),
                                  ))
print summary

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
fig = array_image(rows, cols, ratio, num=1)

# M-A plot --------------------------------------------------------------------
fig = MA_plot(red, green, num=2)

# Profile plot ----------------------------------------------------------------
fig = cgh_profile(locus[locus_order]['start'], ratio[locus_order])

#print xmin, xmax
#plt.axis([xmin, xmax, 0.5, 1.5])
plt.show()



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
