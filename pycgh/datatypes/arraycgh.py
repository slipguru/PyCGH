# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import itertools as it

import numpy as np
from numpy.lib import recfunctions

from ..utils import _file_handle
from cytobands import _chr2int

class ArrayCGH(object):
    

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
        """
        Main data structure representing a generic Array CGH.
        
        
        Parameters
        ----------
        
        id : array_like
            An array containing the name of all probes.
            
        row : array_like
            An array containing, for each probe, the index of the row where the probe lies on the chip.
        
        col : array_like
            An array containing, for each probe, the index of the column where the probe lies on the chip.
        
        reference_signal : array_like
            The intensity values referring to the reference tissue.
        
        test_signal : array_like
            The intensity values referring to the test tissue.
        
        chromosome : array_like
            For each probe, the chromosome where that particular DNA fragment can be found.
        
        start_base : array_like
            The index of the first base of the DNA fragment of each probe.
        
        end_base : array_like
            The index of the last base of the DNA fragment of each probe.
        
        mask : array_like
            An array which indicates for each probe if for some reason did not pass a quality test or is not a real probe (e.g. is a control probe).
           
        """
        
        """
        #. The instantiation of the class needs a reduced number of parameter
        #. The class methods `load` and `save` offer a way to easily import and
           export a processed Array CGH.
        """
        
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
        """
        Loads an aCGH stored in a file and returns a :class:`~pycgh.datatypes.ArrayCGH` object.
        
        Parameters
        ----------
        
        acgh_file : file or str
            Either a file handle or the path to the file containing the aCGH data in text format.
            
        fields : dict, optional (default: ``None``)
            A dictionary which maps names of properties of a :class:`~pycgh.datatypes.ArrayCGH` object (keys) to names of fields found in the input file (values).
        
        delimiter : str, optional (default: ``', '``)
            The string used to separate columns in the CSV file.
        
        missing_values : str, optional (default: ``'--'``)
            The set of strings corresponding to missing data.
        
        Returns
        -------
        
        aCGH : :class:`pycgh.datatypes.ArrayCGH`
            The object representing the Array CGH loaded from the input file.
        """
        
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
        """
        Load an aCGH from a file previously saved using :meth:`pycgh.datatypes.ArrayCGH.save`
        
        Parameters
        ----------
        
        acgh_file : file or str
        
        fields : dict, optional (default: ``None``)
            A dictionary which maps names of properties of a :class:`~pycgh.datatypes.ArrayCGH` object (keys) to names of fields found in the input file (values).
            
        Returns
        -------
        
        aCGH : :class:`pycgh.datatypes.ArrayCGH`
            The object representing the Array CGH.
        """
        
        data = np.load(acgh_file)

        # Try compressed file (default)
        try:
            data = data['acgh']
        except ValueError:
            pass

        mandatory, optionals = _load(data, fields)
        return ArrayCGH(*mandatory, **optionals)

    def savetxt(self, path, fmt="%s", delimiter=', '):
        """
        Saves the :class:`~pycgh.datatypes.ArrayCGH` object as a CSV file using :func:`numpy.savetxt`.
        
        Parameters
        ----------
        
        path : str
            The path of the file which will be saved.
        
        fmt : str, optional (default: ``"%s"``)
            A single format (e.g. ``'%10.5f'``), a sequence of formats, or a multi-format string, e.g. ``'Iteration %d - %10.5f'``, in which case *delimiter* is ignored.
        
        delimiter : str, optional (default: ``', '``)
            The string used to separate columns in the CSV file.
        """
        
        fh = _file_handle(path, mode='w')
        fh.write('%s\n' % delimiter.join(self.names))
        np.savetxt(fh, self._rdata, fmt=fmt, delimiter=delimiter)

    def save(self, path, compressed=True):
        """
        Saves the :class:`~pycgh.datatypes.ArrayCGH` object in the (optionally compressed) numpy binary format (``.npy``/``.npz``).
        
        Parameters
        ----------
        
        path : str
            The path of the file which will be saved.
        
        compressed : bool, optional (default: ``True``)
            Whether to save the file in a compressed (`.npz`) format or in the binary :mod:`numpy` format (``.npy``)
        """
        
        fh = _file_handle(path, mode='w')
        if compressed:
            #np.savez_compressed(fh, acgh=self._rdata)
            np.savez(fh, acgh=self._rdata)
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
