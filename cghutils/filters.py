import numpy as np

# Probes filtering ------------------------------------------------------------
def probes_filter(acgh, filter_controls=True, filter_unmapped=True):

    size = acgh.stat('TotalNumFeatures')
    controls = unmapped = np.zeros(size, dtype=bool) # no filtering

    if filter_controls:
        controls = np.logical_not(acgh.feature('ControlType') == 0)

    if filter_unmapped:
        unmapped = (acgh.feature('SystematicName') == 'unmapped')

    # meaning -> NOT (controls OR unmapped)
    return np.logical_not(np.logical_or(controls, unmapped))


# Conversions algoritm --------------------------------------------------------
def _split_locus(location):
    chr, interval = location.split(':')
    start, end = (int(x) for x in interval.split('-'))

    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr.isdigit():
        chr = chr.zfill(2)

    return chr, start, end

def split_locus_mapping(locations):
    loc_tuples = [_split_locus(x) for x in locations]
    rec_type = np.dtype([('chromosome', 'U20'),
                         ('start', np.int),
                         ('end', np.int),
                        ])
    return np.array(loc_tuples, dtype=rec_type)
