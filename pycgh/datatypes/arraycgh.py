# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import itertools as it

import numpy as np
from numpy.lib import recfunctions

from ..utils import _file_handle
from cytobands import _chr2int

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

    # Common constant
    MISSING_INT = -1
    MISSING_FLOAT = np.nan
    MISSING_STRING = '--'

    def __init__(self, id, row, col, reference_signal, test_signal,
                 chromosome, start_base, end_base, mask=None, **kwargs):

        # Force type for some arguments
        row = np.asanyarray(row, dtype=np.int)
        col = np.asanyarray(col, dtype=np.int)
        reference_signal = np.asanyarray(reference_signal, dtype=np.float)
        test_signal = np.asanyarray(test_signal, dtype=np.float)
        start_base = np.asanyarray(start_base, dtype=np.int)
        end_base = np.asanyarray(end_base, dtype=np.int)
        # Mask indicates probes to hide... default all False

        # Chromosome conversion
        def _conversion(chr):
            chr = str(chr).strip()
            if chr in [str(x) for x in xrange(1, 25)]:
                return int(chr)
            elif chr == 'X':
                return 23
            elif chr == 'Y':
                return 24
            else:
                return self.MISSING_INT

        chromosome = np.asanyarray([_conversion(chr)
                                    for chr in chromosome],
                                   dtype=np.int)

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

        self._rdata =  np.rec.fromarrays(buffer, names=names).view(np.ndarray)
        self._rdata.sort(order=['chromosome', 'start_base'])

        # filtered and masked proxies
        self.F = _ProxyFilter(self.filtered)
        self.M = _ProxyFilter(self.masked)

    def __len__(self):
        """ Raw lenght """
        return len(self._rdata)

    @property
    def size(self):
        """ Valid lenght (equal to __len__ if compacted) """
        # More efficient than inverting the mask
        return len(self._rdata) - (self._rdata['mask']).sum()

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

    def shrink(self):
        self._rdata = self._rdata[~self._rdata['mask']]

    def __setitem__(self, key, value):
        value = np.asanyarray(value)

        # Filtered data handling (autofilled with 0)
        if len(value) == self.size:
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
    def loadtxt(acgh_file, fields=None, delimiter=',', missing_values='--'):
        fh = _file_handle(acgh_file)
        data = np.genfromtxt(fh, delimiter=delimiter, comments='#',
                             names=True, dtype=None, case_sensitive='lower',
                             autostrip=True, invalid_raise=True,
                             missing_values=missing_values,
                             filling_values={   # Default missing values
                                'id': ArrayCGH.MISSING_STRING,
                                'row': ArrayCGH.MISSING_INT,
                                'col': ArrayCGH.MISSING_INT,
                                'reference_signal': ArrayCGH.MISSING_FLOAT,
                                'test_signal': ArrayCGH.MISSING_FLOAT,
                                'chromosome': ArrayCGH.MISSING_INT,
                                'start_base': ArrayCGH.MISSING_INT,
                                'end_base': ArrayCGH.MISSING_INT}
                             )

        # TODO: tweaking performances! Eg. Loadtxt

        mandatory, optionals = _load(data, fields)
        return ArrayCGH(*mandatory, **optionals)

    @staticmethod
    def load(acgh_file, fields=None):
        fh = _file_handle(acgh_file)
        data = np.load(fh)

        # Try compressed file (default)
        try:
            data = data['acgh']
        except ValueError:
            pass

        mandatory, optionals = _load(data, fields)
        return ArrayCGH(*mandatory, **optionals)

    def savetxt(self, path, fmt="%s", delimiter=', '):
        fh = _file_handle(path, mode='w')
        fh.write('%s\n' % delimiter.join(self.names))
        np.savetxt(fh, self._rdata, fmt=fmt, delimiter=delimiter)

    def save(self, path, compressed=True):
        fh = _file_handle(path, mode='w')
        if compressed:
            np.savez_compressed(fh, acgh=self._rdata)
        else:
            np.save(fh, self._rdata)


def _load(data, fields):
    # keys: expected aCGH fields; values: input fields
    file_fields = dict((v, v) for v in data.dtype.names)
    if not fields is None:
        for name in fields:
            fields[name] = fields[name].lower()
            del file_fields[fields[name]]
        file_fields.update(fields)

    optional_columns = set(file_fields.keys()).difference(ArrayCGH.COL_NAMES)
    optional_inputs = dict((k, data[file_fields[k]]) for k in optional_columns)

    mandatory_inputs = [data[file_fields[k]] for k in ArrayCGH.COL_NAMES]

    return mandatory_inputs, optional_inputs


class _ProxyFilter(object):
    def __init__(self, getter):
        self._getter = getter # reference to main object

    def __getitem__(self, key):
        return self._getter(key)
