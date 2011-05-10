import itertools as it
import csv

import numpy as np
from numpy.lib import recfunctions

class ArrayCGH(object):

    COL_NAMES = ('id', 'row', 'col',
                 'reference_signal', 'test_signal',
                 'chromosome', 'start_base', 'end_base')

    def __init__(self, id, row, col, reference_signal, test_signal,
                 chromosome, start_base, end_base, mask=None, **kwargs):

        # Default: no optional inputs
        buffer = [id, row, col, reference_signal, test_signal,
                  chromosome, start_base, end_base]
        names = list(self.COL_NAMES)

        # Read mask as a standard optional parameter
        if mask is None:
            # std "False-mask"
            kwargs['mask'] = np.zeros(len(buffer[0]), dtype=np.bool)
        else:
            kwargs['mask'] = np.asanyarray(mask)

        # Extend the list of inputs
        if kwargs:
            kwnames, kwvalues = zip(*kwargs.items())
            buffer = it.chain(buffer, kwvalues)
            names.extend(kwnames)

        if len(set(names)) < len(names):
            raise ValueError('optional parameters name duplication')

        self._rdata =  np.rec.fromarrays(buffer, names=names).view(np.ndarray)

    def __len__(self):
        """ mumble mumble """
        return len(self._rdata)

    @property
    def names(self):
        return self._rdata.dtype.names

    def __getitem__(self, key):
        """ Returns a copy of the filtered data """
        return self._rdata[key][~self._rdata['mask']]

    def __setitem__(self, key, value):
        value = np.asanyarray(value)

        # Filtered data handling
        if len(value) == (~self._rdata['mask']).sum():
            full_value = np.zeros(len(self._rdata), dtype=value.dtype)
            full_value[~self._rdata['mask']] = value
            value = full_value
        elif len(value) != len(self._rdata):
            raise ValueError('wrong dimension')

        if key in self.names: #Update
            self._rdata[key] = value
        else: #Add
            self._rdata = recfunctions.append_fields(self._rdata, names=key,
                                                     data=value, usemask=False)

    def unfiltered(self, key, copy=False):
        """ Returns a copy if copy=True """
        if copy:
            return self._rdata[key].copy()
        else:
            return self._rdata[key]

    def masked(self, key, copy=False, fill_value=None):
        """ Returns a copy if copy=True """
        return np.ma.masked_array(self._rdata[key], self._rdata['mask'],
                                  copy=copy, fill_value=fill_value)

    def sort(self, order):
        self._rdata.sort(order=order)

    def __repr__(self):
        return repr(self._rdata)

    def __str__(self):
        return str(self._rdata)

    @staticmethod
    def load(path, fields=None, delimiter=','):
        data = np.genfromtxt(path, delimiter=delimiter,
                             names=True, dtype=None, case_sensitive='lower',
                             autostrip=True, invalid_raise=True)

        # chiavi sono campi aCHG, valori cosa leggere da input
        file_fields = dict((v, v) for v in data.dtype.names)
        if not fields is None:
            for name in fields:
                fields[name] = fields[name].lower()
                del file_fields[fields[name]]
            file_fields.update(fields)

        optional_columns = set(file_fields.keys()).difference(ArrayCGH.COL_NAMES)
        optional_inputs = dict((k, data[file_fields[k]]) for k in optional_columns)

        return ArrayCGH(*[data[file_fields[k]] for k in ArrayCGH.COL_NAMES],
                        **optional_inputs)

    def save(self, path, delimiter=','):
        with open(path, 'wb') as out:
            writer = csv.writer(open(path, 'wb'), delimiter=delimiter)
            writer.writerow(self.names)
            for a in self._rdata:
                writer.writerow(a)
