import itertools as it
import csv

import numpy as np

class ArrayCGH(object):

    COL_NAMES = ('id', 'row', 'col',
                 'reference_signal', 'test_signal',
                 'chromosome', 'start_base', 'end_base')

    def __init__(self, input, mask=None, **kwargs):

        # Default: no optional inputs
        buffer = list(input)
        names = list(self.COL_NAMES)

        if len(names) != len(buffer):
            raise ValueError('missing mandatory inputs, only %d provided' % len(buffer))

        # Read mask as a standard optional parameter
        if not mask is None:
            kwargs['mask'] = mask

        # Extend the list of inputs
        if kwargs:
            kwnames, kwvalues = zip(*kwargs.items())
            buffer = it.chain(buffer, kwvalues)
            names.extend(kwnames)

        if len(set(names)) < len(names):
            raise ValueError('optional parameters name duplication')

        self._rdata =  np.rec.fromarrays(buffer, names=names).view(np.ndarray)

    def __len__(self):
        return len(self._rdata)

    @property
    def names(self):
        return self._rdata.dtype.names

    def __getitem__(self, key):
        return self._rdata[key]

    def __repr__(self):
        return repr(self._rdata)

    def __str__(self):
        return str(self._rdata)

    @staticmethod
    def load(path, fields=None, delimiter=','):
        data = np.genfromtxt(path, delimiter=delimiter,
                             names=True, dtype=None, case_sensitive='lower',
                             autostrip=True, invalid_raise=True)

        # chiavi sono campi aCHG, valori cosa leggere da tmp
        file_fields = dict((v, v) for v in data.dtype.names)
        if not fields is None:
            for name in fields:
                fields[name] = fields[name].lower()
                del file_fields[fields[name]]
            file_fields.update(fields)

        optional_columns = set(file_fields.keys()).difference(ArrayCGH.COL_NAMES)
        optional_inputs = dict((k, data[file_fields[k]]) for k in optional_columns)

        return ArrayCGH([data[file_fields[k]] for k in ArrayCGH.COL_NAMES],
                        **optional_inputs)

    def save(self, path, delimiter=','):
        with open(path, 'wb') as out:
            writer = csv.writer(open(path, 'wb'), delimiter=delimiter)
            writer.writerow(self.names)
            for a in self._rdata:
                writer.writerow(a)
