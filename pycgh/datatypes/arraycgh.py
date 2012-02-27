# -*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import itertools as it

import numpy as np
from numpy.lib import recfunctions

from ..io.utils import _file_handle

class ArrayCGH(object):
    """ Main data structure repersenting a generic Array CGH.

    Features
    --------
    1. The instantiation of the class needs a reduced number of parameter
    2. The class methods `load` and `save` offer a way to easily import and
       export a processed Array CGH.

    """

    # Mandatory column names
    COL_NAMES = ('id', 'row', 'col',
                 'reference_signal', 'test_signal',
                 'chromosome', 'start_base', 'end_base')
    # 'mask' has a special meaning but it is not mandatory

    def __init__(self, id, row, col, reference_signal, test_signal,
                 chromosome, start_base, end_base, mask=None, **kwargs):

        # Mask indicates probes to hide... default all False

        # Default: no optional inputs
        buffer = [id, row, col, reference_signal, test_signal,
                  chromosome, start_base, end_base]
        names = list(self.COL_NAMES)

        # Read mask as a standard optional parameter
        if mask is None:
            # std "False-mask"
            kwargs['mask'] = np.zeros(len(buffer[0]), dtype=np.bool)
        else:
            kwargs['mask'] = np.asanyarray(mask, dtype=np.bool)

        # Extend the list of inputs
        if kwargs:
            kwnames, kwvalues = zip(*kwargs.items())
            buffer = it.chain(buffer, kwvalues)
            names.extend(kwnames)

        if len(set(names)) < len(names):
            raise ValueError('optional parameters name duplication')

        self._rdata =  np.rec.fromarrays(buffer, names=names).view(np.ndarray)

        # We have to evaluate the efficency to have a record array
        # or to have a dictionary of column vector...
        # With a record array we can also access a Probe row... is it useful?
        # In order to have this feature, are we sacrifying efficiency
        # accessing columns??

        # Access by probe is probably not so useful... than the rec array,
        # could be a simple interface to with operate

    def __len__(self):
        """ mumble mumble """
        return len(self._rdata)

    @property
    def names(self):
        return self._rdata.dtype.names

    # - raw data (view) - Default
    # - masked data (view or copy)
    # - filtered data (copy)
    # Default only as a copy?????

    def __getitem__(self, key):
        """ Returns a view of the raw data """
        return self._rdata[key].view()

    def masked(self, key, copy=False, fill_value=None):
        """ Returns a copy if copy=True """
        return np.ma.masked_array(self._rdata[key], self._rdata['mask'],
                                  copy=copy, fill_value=fill_value)

    def filtered(self, key):
        """ Returns a copy """
        return self._rdata[~self._rdata['mask']][key]

    def sort(self, order):
        self._rdata.sort(order=order)

    def __setitem__(self, key, value):
        value = np.asanyarray(value)

        # Filtered data handling (autofilled with 0)
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

    def __repr__(self):
        return repr(self._rdata)

    def __str__(self):
        return str(self._rdata)

    @staticmethod
    def load(acgh_file, fields=None, delimiter=','):
        fh = _file_handle(acgh_file)
        data = np.genfromtxt(fh, delimiter=delimiter, comments='#',
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
        fh = _file_handle(path, mode='w')
        fh.write('%s\n' % delimiter.join(self.names))
        np.savetxt(fh, self._rdata, fmt="%s", delimiter=delimiter)

