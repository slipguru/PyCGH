#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

from pycgh import PyCGHException
from ..utils import _file_handle


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

    #The class encapsulates the concepts of a collection of omogeneous
    #labeled samples, described by a fixed number of named variables.
    #
    #Another way to handle the class is to imagine it as an abstract
    #representation of a two-entry table, where on the rows we have the
    #different labeled sample and on the columns the different named variables.
    #
    #In the context of the class, the the two terms `label` and `name` are
    #respectively associated with `samples` and `variables`
    #(eg: "the sample label" and "the variable name").
    #
    #Moreover, each variables has also a specific type. The internal
    #representation of the data inherits the data type representation from
    #the numpy.ndarray data structure.

class DataTable(object):
    def __init__(self, data, rlabels=None, clabels=None, dtype=None):
        """ Data are not copied if already an ndarray! """

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
            return [str(x).strip() for x in labels]

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
        """
        In input we can have:
            - A Slice
            - A Label (string-like)
            - An integer (index)
            - An iterable of labels or integers

            - A tuple (pair) of the above types
        """

        if isinstance(item, tuple):
            if len(item) == 1: # special case
                ritem = self._map_item(item[0], self._rdict)
                citem = None
            elif len(item) == 2:
                ritem = self._map_item(item[0], self._rdict)
                citem = self._map_item(item[1], self._cdict)
            else:
                raise PyCGHException('maximum 2 dimension allowed')
        else:
            ritem = self._map_item(item, self._rdict)
            citem = None

        try:
            if citem is None:
                return self._data[ritem]
            return self._data[ritem, citem]
        except (ValueError, IndexError):
            raise PyCGHException('invalid labels names or '
                                 'indexes: %s' % str(item))


    def _map_item(self, item, d): # Exception postoponed
        if isinstance(item, slice):
           return slice(d.get(item.start, item.start),
                        d.get(item.stop,  item.stop),
                        d.get(item.step,  item.step))
        elif str(item) == item: # string-like
            return d.get(item, item)
        else:
            try: # iterable (eventually mixed)
                return [d.get(i, i) for i in item]
            except TypeError:
                return d.get(item, item) # index

    def rindex(self, label):
        try:
            return self._rdict[label]
        except KeyError:
            raise PyCGHException('row label not found: %s' % label)

    def cindex(self, label):
        try:
            return self._cdict[label]
        except KeyError:
            raise PyCGHException('column label not found: %s' % label)

    def __array__(self):
        """ Ensures data are not copied """
        return self._data

    @property
    def rnum(self):
        return len(self._rdict)

    @property
    def cnum(self):
        return len(self._cdict)

    def __len__(self):
        return self.rnum

    @property
    def rlabels(self):
        return self._rlabels

    @property
    def clabels(self):
        return self._clabels

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def shape(self):
        return self._data.shape

    @staticmethod
    def load(dt_file, delimiter='\t', dtype=float):
        fh = _file_handle(dt_file)

        try:
            data = np.loadtxt(fh, delimiter=delimiter, dtype='S')
            clabels = data[0,1:]
            rlabels = data[1:,0]
        except:
            raise PyCGHException('failed loading file %s' % dt_file)

        return DataTable(data[1:, 1:], rlabels, clabels, dtype)

    def save(self, path, fmt="%s", delimiter=', '):
        fh = _file_handle(path, mode='w')
        fh.write('%s%s' % ('*', delimiter))
        fh.write('%s\n' % delimiter.join(self.clabels))
        formatter = '%s\n' % delimiter.join((fmt,) * len(self.clabels))

        for line, label in zip(self._data, self._rlabels):
            fh.write('%s%s' % (label, delimiter))
            fh.write(formatter % tuple(line))
