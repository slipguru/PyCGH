from ..utils import _file_handle, location_normalization
from ..datatypes.arraycgh import ArrayCGH

# filter valid return only valid probes
def ucsc_mapping(ucsc_file, filter_valid=False):
    ucsc_mapping = dict()
    for line in _file_handle(ucsc_file):
        _, chr, sb, eb, id = line.split()
        value = location_normalization(chr, int(sb)+1, int(eb))
        if not filter_valid:
            ucsc_mapping[id] = value
        elif not value[3]: # Only If valid
            ucsc_mapping[id] = value[:3]

    return ucsc_mapping

