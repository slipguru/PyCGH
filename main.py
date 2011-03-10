from collections import defaultdict

import numpy as np
from matplotlib import pylab as plt
from matplotlib import mlab

from cghutils.readers import AgilentReader
from cghutils.filters import probes_filter, split_locus_mapping

# GPL5477
ABS_PATH = '/home/sabba/Phd/Tonini_IST/Group1/1115_G1_251495014936_1_2.txt'

acgh = AgilentReader(ABS_PATH)

# Remove controls, unmapped and strange chr and get locations
valid_probes = probes_filter(acgh)
locations = acgh.feature('SystematicName')[valid_probes]

#ratio = acgh.feature('LogRatio')[valid_probes]
#green = acgh.feature('gMedianSignal')[valid_probes]
#red = acgh.feature('rMedianSignal')[valid_probes]

# Extact locus and sort by them
locus = split_locus_mapping(locations)
locus_order = np.argsort(locus)

# Print a chromosome summary
summary = mlab.rec_groupby(locus,
                           groupby=('chromosome',),
                           stats=(('start', np.min, 'from'),
                                   ('end', np.max, 'to'),
                                  ))
print summary






#for item in summary:
#    if item['clones'] > 1:
#        print item



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

## Collect duplication
#counter = defaultdict(list)
#for s, l in zip(signals, locations):
#    counter[l].append(s)
#
## Plot Boxplots of duplicated probes
#duplicated = list()
#for c in counter:
#    if len(counter[c]) > 1:
#        duplicated.append(counter[c])
#
#plt.boxplot(duplicated)
