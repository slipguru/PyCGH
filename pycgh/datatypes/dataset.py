#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

class Dataset(object):
    """Dataset data type.
    
    The class encapsulate the concepts of a collection of omogeneous
    labeled samples, described by a fixed number of named variables.
    
    Another way to handle the class is to imagine it as an abstract
    representation of a two-entry table, where on the rows we have the
    different labeled sample and on the columns the different named variables.
    
    In the context of the class, the the two terms `label` and `name` are
    respectively associated with `samples` and `variables`
    (eg: "the sample label" and "the variable name").
    
    Moreover, each variables has also a specific type. The internal
    representation of the data inherits the data type representation from
    the numpy.ndarray data structure.
    
    """
    
    def __init__(self, names, types=np.float32,
                 data=None, labels=None):
        """ Dataset constructor.
        
        The class constructor requires only the description of the samples
        as list of variable `names` and `types`.
        The input variable `types` may be a simple python (or numpy) type
        and the class assumes that all the variables have the same type
        (as default types is `np.float32`).
        
        Data has to be a list of samples (in a form convertible in a
        numpy 2d ndarray, eg. list of lists, list of 1d ndarrays ...),
        where the lenght of each sample must be equal the lenght of `names`.
        
        Finally, `labels' is the list of sample labels. If it is None,
        the class automatically will use the labels ['sample0', 'sample1', ...].
        """
        
        ##TODO: if types is a list we have to build consistently _rdata
        
        self._names = dict((n, i) for i, n in enumerate(names))

        if data is None:
            self._rdata = np.empty((0, len(names)), dtype=types)
            self._samples = {}
        else:
            data = np.asanyarray(data, dtype=types)
            r, c = data.shape
            if len(names) != c:
                raise ValueError('data column number must be equal to number of column names')

            self._rdata = data
            if labels is None:
                self._samples = dict(('sample%d' % i, i) for i in xrange(1, r+1))
            else:
                self._samples = dict((name, i) for i, name in enumerate(labels))

    def index_of(self, key):
        return self._samples[key]

    ###############TO TEST
    def asarray(self):
        return np.asanyarray(self._rdata)

    @property
    def samples(self):
        return self._samples.keys()

    @property
    def names(self):
        return self._names.keys()

    def append(self, key, values):
        " ADD "
        values = np.asanyarray(values)
        values.shape = (1, len(values))

        if key is None:
            key = 'sample%d' % (max([-1] + self._samples.values()) + 1)

        if key in self._samples:
            raise ValueError('sample already present')

        self._rdata = np.r_[self._rdata, values]
        self._samples[key] = (len(self._rdata) - 1)

    def __getitem__(self, key):
        " READ "
        if key in self._samples:
            return self._rdata[self._samples[key]]
        return self._rdata[key]

    def __setitem__(self, key, value):
        " CHANGE "
        if not key in self._samples:
            raise ValueError('sample not existent')

        values = np.asanyarray(value)
        values.shape = (1, len(values))
        self._rdata[self._samples[key]] = values

    def __delitem__(self, key):
        " DELETE not efficient "
        if not key in self._samples:
            raise ValueError('sample not existent')

        # Get the sample name of the last row
        from operator import itemgetter as ig
        last_key, last_index = max(self._samples.iteritems(), key=ig(1))

        # Copy the last row in the row to eliminate and delete the last one
        rm_index = self._samples[key]
        self._rdata[rm_index] = self._rdata[last_index]
        self._rdata = np.delete(self._rdata, last_index, 0)

        # Update the mapping between sample names and indexes
        self._samples[last_key] = rm_index
        del self._samples[key]

    def __len__(self):
        return len(self._rdata)

    def __str__(self):
        return str(self._rdata)

    @staticmethod
    def load(file_path, delimiter='\t'):
        # For some reason genfromtxt removes the '.' in the names string
        f = open(file_path)
        header = '#'
        while header.startswith('#'):
            header = f.readline()
        names = header.strip().split(delimiter)
        attrs = names[1:]

        data = np.genfromtxt(f, dtype=None,
                             delimiter=delimiter, names=names)
        dnames = data.dtype.names

        samples = data[names[0]]
        data = (data[list(dnames[1:])].view(dtype=np.float)
                                      .reshape((len(samples), len(attrs)))
               )

        return LabeledMatrix(names=attrs,
                             data=data,
                             labels=samples)

    def save(self, path, delimiter='\t'):
        from operator import itemgetter as ig
        import csv

        with open(path, 'wb') as out:
            writer = csv.writer(open(path, 'wb'), delimiter=delimiter)

            names = [n for n,i in sorted(self._names.iteritems(), key=ig(1))]
            writer.writerow(['sample'] + names)

            for s in self._samples:
                idx = self._samples[s]
                writer.writerow([s] + self._rdata[idx].tolist())