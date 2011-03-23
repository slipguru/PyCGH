import numpy as np

# Probes filtering ------------------------------------------------------------
def probes_filter(acgh, filter_controls=True, filter_unmapped=True):

    try:
        size = acgh.stat('TotalNumFeatures')
    except:
        size = acgh.feature('SystematicName').size
    controls = unmapped = np.zeros(size, dtype=bool) # no filtering

    if filter_controls:
        controls = np.logical_not(acgh.feature('ControlType') == 0)

    if filter_unmapped:
        unmapped = (acgh.feature('SystematicName') == 'unmapped')

    # meaning -> NOT (controls OR unmapped)
    return np.logical_not(np.logical_or(controls, unmapped))


# Conversions algoritm --------------------------------------------------------
def _split_mapping(location):
    chr, interval = location.split(':')
    start, end = (int(x) for x in interval.split('-'))

    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr.isdigit():
        chr = chr.zfill(2)

    return chr, start, end

def split_mappings(locations):
    """ Converts a list of locations strings in a numpy record array
    indexed by 'chromosome', 'start', 'end' """
    loc_tuples = [_split_mapping(x) for x in locations]
    rec_type = np.dtype([('chromosome', 'U2'),
                         ('start', np.int),
                         ('end', np.int),
                        ])
    return np.array(loc_tuples, dtype=rec_type)
