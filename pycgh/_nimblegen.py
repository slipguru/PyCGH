#-*- coding: utf-8 -*-
import numpy as np

from _arrayCGH import ArrayCGH


TYPE_MAP = {'IMAGE_ID': unicode,
            'GENE_EXPR_OPTION': unicode,
            'SEQ_ID': unicode,
            'PROBE_ID': unicode,
            'POSITION': int,
            'X': int,
            'Y': int,
            'MATCH_INDEX': int,
            'SEQ_URL': unicode,
            'PM': float,
            'MM': float}

def _read_all_data(fileobj, header, delimiter):
    out = dict()
    for line in fileobj:
        data = line.strip().split(delimiter)
        for h, d in zip(header, data):
            out.setdefault(h, []).append(TYPE_MAP[h](d))

    # Conversion in numpy array
    for k in out:
        out[k] = np.asarray(out[k], dtype=TYPE_MAP[k])

    return out

def _read_pair(file_path, delimiter):
    pair_file = open(file_path, 'r')
    meta_line = pair_file.readline()[1:].strip().split(delimiter)
    meta = dict((k, v) for (k, v) in (x.split('=') for x in meta_line))

    header = pair_file.readline().strip().split(delimiter)
    data = _read_all_data(pair_file, header, delimiter)

    return meta, data

def _split_mapping(location):
    if 'CHR' in location:
        chr, start = location.split('00P')
        start = int(start)
        try:
            chr = int(chr[3:])
        except ValueError:
            if 'X' in chr:
                chr = 23
            elif 'Y' in chr:
                chr = 24
            else:
                raise ValueError('unrecognized chromosome position '
                                 '(probe id: %s)' % location)
        return (chr, start, False)
    else:
        return (NimblegenCGH.INVALID_INT, NimblegenCGH.INVALID_INT, True)

class NimblegenCGH(ArrayCGH):

    INVALID_INT = -9999

    def __init__(self, *args, **kwargs):
        return super(NimblegenCGH, self).__init__(*args, **kwargs)

    @staticmethod
    def load(test_path, reference_path, delimiter='\t'):
        test_meta, test_data = _read_pair(test_path, delimiter)
        ref_meta, ref_data = _read_pair(reference_path, delimiter)

        if not np.all(test_data['PROBE_ID'] == ref_data['PROBE_ID']):
            raise ValueError('probe id mismatch between pair files')

        locations = zip(*(_split_mapping(x) for x in test_data['PROBE_ID']))
        chromosome, start_base, mask = locations

        # Creation and dinamyc attachment of useful informations
        aCGH = NimblegenCGH(id=test_data['PROBE_ID'],
                            row=test_data['Y'],
                            col=test_data['X'],
                            reference_signal=ref_data['PM'],
                            test_signal=test_data['PM'],
                            chromosome=chromosome,
                            start_base=start_base,
                            end_base=start_base, #TODO,
                            mask=mask)
        aCGH.META_TEST = test_meta
        aCGH.META_REF = ref_meta
        aCGH.NAMES_MAP = dict(zip(ArrayCGH.COL_NAMES[:6],
                                  ('PROBE_ID', 'Y', 'X', 'PM', 'PM', 'PROBE_ID'))
                             )

        return aCGH



