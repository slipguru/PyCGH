from ..utils import _file_handle
from ..datatypes.arraycgh import ArrayCGH

# filter valid return only valid probes
def ucsc_mapping(ucsc_file, filter_valid=False):
    ucsc_mapping = dict()
    for line in _file_handle(ucsc_file):
        _, chr, sb, eb, id = line.split()
        value = _split_mapping(chr, int(sb)+1, int(eb))
        if not filter_valid:
            ucsc_mapping[id] = value
        elif not value[3]: # Only If valid
            ucsc_mapping[id] = value[:3]

    return ucsc_mapping

def _split_mapping(chr, start, end):
    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr == 'X':
        return 23, start, end, False
    elif chr=='Y':
        return 24, start, end, False
    else:
        try:
            return int(chr), start, end, False
        except ValueError: # unplaceable probe eg chrUn
            return (ArrayCGH.MISSING_INT,
                    ArrayCGH.MISSING_INT,
                    ArrayCGH.MISSING_INT, True)
