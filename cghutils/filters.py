import numpy as np

# Probes filtering ------------------------------------------------------------
def probes_filter(acgh, filter_controls=True, filter_unmapped=True,
                        filter_strange_chr=True):

    size = acgh.stat('TotalNumFeatures')
    controls = unmapped = strange_chr = np.zeros(size, dtype=bool) # no filtering

    if filter_controls:
        controls = np.logical_not(acgh.feature('ControlType') == 0)

    if filter_unmapped:
        unmapped = (acgh.feature('SystematicName') == 'unmapped')


    # Probabilmente possiamo mantenere questi cromosomi eliminando
    # il suffisso. Questo si riferisce al mapping della sequenza con il
    # relativo gene che nelle nostre analisi p un problema secondario!
    if filter_strange_chr: #e.g. _random
        strange_chr = (np.char.find(acgh.feature('SystematicName'), '_') != -1)

    # meaning -> NOT (controls OR unmapped OR strange_cgh)
    return np.logical_not(np.logical_or(np.logical_or(controls, unmapped),
                                        strange_chr))


# Conversions algoritm --------------------------------------------------------
def _split_locus(location):
    chr, interval = location.split(':')
    start, end = (int(x) for x in interval.split('-'))

    chr = chr.replace('chr', '')
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
