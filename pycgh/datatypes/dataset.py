#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np
from pycgh import PyCGHException

### L'utilità di questa classe dovrebbe essere limitata alla
### possibilità sia di leggere tabelle cliniche (anche se basta il csv reader,
### se non fosse che è un reader e non un parser che ti mette tutto in memoria.)
### rispetto al csv reader permette slicing su righe e colonne offreddo accesso
### alle righe, quindi ha senso di esistere.
###
### Un secondo uso è quello di gestire dataset si sample viste sempre come
### tabelle csv e offrendo le stesse funzionalità
### oltre che alla possibilità di aggiungere ed eliminare colonne e righe
### e risalvare il file tutto molto più semplice rispetto ad un csv reader.
###
### Sul tipo, per ora limitiamoci ai tipi base, non ai record array.
### Valgono solo i ".isbuiltin == 1", così posso salvarla come csv
###
### Quindi l'astrazione corretta è quella di Table con row and columns

### Semplificare avendo internamente due indici (mappe) sia per colonne
### che per righe. Accetto in input un dtype numpy...

### Di cosa ho bisogno: lettura scrittura csv, non dimenticarmi delle etichette
### larry sarebbe più che sufficiente se non fosse che legge una merda i csv.

class DataTable(object):
    def __init__(self, data, rlabels=None, clabels=None, dtype=None):

        ### TODO: check data dimension == 2
        ### Pur rlabels and clabels into 2 dict

        self._data = self._create_data(data, dtype)
        r, c = self._data.shape
        self._rlabels = self._autolabels(rlabels, r, 'r')
        self._clabels = self._autolabels(clabels, c, 'c')

        self._rdict = dict(zip(self._rlabels, xrange(len(self._rlabels))))
        self._cdict = dict(zip(self._clabels, xrange(len(self._clabels))))

    def _create_data(self, data, dtype):
        # Because, if None the type will be automatically inferred by numpy
        if dtype is None:
            return np.asanyarray(data)
        else:
            dtype = np.dtype(dtype)
            if dtype.isbuiltin:
                return np.asanyarray(data, dtype=dtype)
            else:
                raise PyCGHException('only builtin dtype permitted.')

    def _autolabels(self, labels, dim, prefix):
        if labels is None:
            return ['%s%d' % (prefix, i) for i in xrange(dim)]

        if len(labels) == dim:
            return [str(x) for x in labels]

        raise PyCGHException('wrong %slabels dimension.' % prefix)

    def __str__(self):
        return '\n'.join(('rows: %s' % ', '.join(self.rlabels),
                          'cols: %s' % ', '.join(self.clabels),
                          'data:\n%s' % str(self._data)))

    def _get_indexes(self, item):
        # Row labels access
        if str(item) == item:
            return item, None

        # Index or Tuple access
        try:
            ritem, citem = item
        except TypeError:
            # Not Iterable object
            ritem, citem = item, None
        except ValueError:
            # Not valid tuple
            if len(item) > 2:
                raise PyCGHException('maximum 2 dimension allowed')
            else:
                # One-element tuple
                ritem, citem = item[0], None

        return ritem, citem

    def __getitem__(self, item):
        # Indexes extraction
        ritem, citem = self._get_indexes(item)

        print 'SABBBBBBBBBBBBBB'

        print
        print item, ":", ritem, citem, '***',

        # Maps or itself
        def map_item(item, d):
            if isinstance(item, slice):
               return slice(d.get(item.start, item.start),
                            d.get(item.stop,  item.stop),
                            d.get(item.step,  item.step))
            elif str(item) == item: # single string
                return d.get(item, item)
            elif isinstance(item, tuple): # iterable (eventually mixed)
                return [d.get(i, i) for i in item]
            else:
                return d.get(item, item)     # index

        ritem = map_item(ritem, self._rdict)
        citem = map_item(citem, self._cdict)

        print ritem, citem

        if citem is None:
            return self._data[ritem]
        return self._data[ritem, citem]

    @property
    def row_num(self):
        return 2

    @property
    def col_num(self):
        return 3

    def __len__(self):
        return self.row_num

    @property
    def rlabels(self):
        return self._rlabels

    @property
    def clabels(self):
        return self._clabels

    @property
    def dtype(self):
        return self._data.dtype

class Dataset(object):
    """Dataset data type.

    The class encapsulates the concepts of a collection of omogeneous
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

    def __init__(self, names, types=float,
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

        Finally, `labels` is the list of sample labels. If it is None,
        the class automatically will use the labels ['sample0', 'sample1', ...].

        TODO: Explain that some types need explicitly the dimension.
        eg. if you want samples with an int, a float and a string of maximum
        10 characters you have to pass a type in numpy format:
        [int, float, 'S10']

        """

        # 'types' may be a list or a single data type
        try:
            if str(types) == types:         # is a single string
                self._types = (types,)*len(names)
            elif len(types) == len(names):  # is a list of types
                self._types = tuple(types)
            else:
                raise PyCGHException("lenght mismatch between "
                                     "'names' and 'types'")
        except TypeError: # if is a single type 'len' raise a TypeError
            self._types = (types,)*len(names)

        # Datatype for the iternal representation of the data (raw np.ndarray)
        datatype = np.dtype({'names': names,
                             'formats': self._types})

        # Raw data creation
        if data is None:
            self._rdata = np.empty((0,), dtype=datatype)
            self._samples = {}
        else:
            try:
                self._rdata = np.asanyarray([tuple(row) for row in data], dtype=datatype)
                if labels is None:
                    self._samples = dict(('sample%d' % i, i)
                                         for i in xrange(1, len(self._rdata)+1))
                elif len(labels) == len(self._rdata):
                    self._samples = dict((name, i) for i, name in enumerate(labels))
                else:
                    raise PyCGHException('data not compatible with the types')
            except ValueError:
                raise PyCGHException("lenght mismatch between "
                                     "'data' and 'labels'")

    def __str__(self):
        return str(self._rdata)

    def __repr__(self):
        return repr(self._rdata)

    # Variables informations --------------------------------------------------
    @property
    def dim(self):
        """Number of variables."""
        return len(self._rdata.dtype)

    @property
    def names(self):
        """Ordered tuple of variables names."""
        return self._rdata.dtype.names

    @property
    def types(self):
        """Ordered tuple of variables types."""
        return self._types

    # Samples informations ----------------------------------------------------
    def __len__(self):
        """Number of samples."""
        return len(self._rdata)

    @property
    def labels(self):
        """Ordered tuple of sample labels."""
        import operator as op

        if self._samples:
            sorted_items = sorted(self._samples.iteritems(),
                                  key=op.itemgetter(1))
            return zip(*sorted_items)[0]
        else:
            return tuple()

    # Samples managing --------------------------------------------------------
    def add(self, label, values):
        """Add the given sample in the dataset collection.

        The position in the matrix is selected for efficiency,
        you can get it by the `index_of` method.
        """
        values = np.asanyarray(tuple(values), dtype=self._rdata.dtype)
        values.shape = (1,)

        if label is None:
            label = 'sample%d' % (len(self._rdata) + 1)

        if label in self._samples:
            raise ValueError('sample already present')

        self._rdata = np.r_[self._rdata, values]
        self._samples[label] = (len(self._rdata) - 1)

    def __getitem__(self, key):
        """Get a sample.

        The user can pass a label or an index.
        """
        # is a label
        try:
            if key in self._samples:
                return self._rdata[self._samples[key]]
        except TypeError:
            pass

        # is an index
        return self._rdata[key]

    def index_of(self, label):
        """Return the index of a specified labeled sample."""
        return self._samples[label]



    ###############TO TEST
    def toarray(self):
        return np.asarray(self._rdata.tolist())






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

        return Dataset(names=attrs,
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
