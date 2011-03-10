import os
from csv import DictReader

import numpy as np

from cghutils.readers import AgilentReader, GPLReader
from cghutils.filters import probes_filter, split_locus_mapping

DATA_PATH = '/home/sabba/Phd/Tonini_IST/'
CGH_PATH = os.path.join(DATA_PATH, 'Group1')
PLATFORM_PATH = os.path.join(DATA_PATH, 'Platforms')

# Reading Clinical data
info_file = open(os.path.join(DATA_PATH, 'info_Row_aCGH.csv'), 'rb')
samples_info = DictReader(info_file, dialect='excel')

for sample in samples_info:

    if sample['Sample name'] == 'Sample 21':
        continue

    platform = sample['platform']
    print
    print '*'*80
    print '%s (%s)...' % (sample['Sample name'], platform)

    # Reading agilent result --------------------------------------------------
    aCGH_path = os.path.join(CGH_PATH,
                             sample['Agilent Feature Extraction file'])
    acgh = AgilentReader(aCGH_path)
    feature_num = acgh.feature('FeatureNum')
    probe_name = acgh.feature('ProbeName')
    systematic_name = acgh.feature('SystematicName')

    # Reading platform informations -------------------------------------------
    platform_path = os.path.join(PLATFORM_PATH, '%s.txt' % platform)
    platform = GPLReader(platform_path)
    id = platform.field('ID')
    spots_id = platform.field('SPOT_ID')
    locations = platform.field('CHROMOSOMAL_LOCATION')

    # to simulate agilent
    for i, (spot, loc) in enumerate(zip(spots_id, locations)):
        if loc == u'':
            locations[i] = spots_id[i]

    agilent_dict = dict(zip(feature_num, zip(systematic_name, probe_name)))
    platform_dict = dict(zip(id, zip(locations, spots_id)))

    assert np.max(feature_num) == np.max(id)

    print
    print '\t'.join(('Feat'.ljust(4),
                     'Locus'.ljust(25),
                     'Probe'.ljust(12),
                     'PlatLocus'.ljust(25),
                     'PlatProbe'.ljust(12)))
    print '-' * 100

    # patform_dict contains some feature to ignore not present
    # on the agilent file
    for id in agilent_dict:
        mapping, probe = agilent_dict[id]
        plat_mapping, plat_probe = platform_dict[id]

        # Same probe but differet mapping
        if probe == plat_probe and mapping != plat_mapping:
            print '\t'.join((str(id).ljust(4),
                             mapping.ljust(25),
                             probe.ljust(12),
                             plat_mapping.ljust(25),
                             plat_probe.ljust(12)))
