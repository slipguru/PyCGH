#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np
from pycgh import PyCGHException

### Potrei implementare la classe in modo che sia
### mono tipo
### In questo modo avrei un classico ndarray
### Questo la differenzierebbe dall'arraycgh anche se... solo l'id è stringa
### .... potrei semplicemente litarmi ai tipi numerici,
### pero se li mischio il problema è o stesso.

### sia larry che pandas sono orientati a serie, non tanto al concetto
### di sample e comunque usano semplici array numpy...

### Semplificare avendo internamente due indici (mappe) sia per colonne
### che per righe. Accetto in input un dtype numpy... potrei ereditare
### La classe funziona come un normalissimo array numpy senza perdere nulla
### della sua efficenza, potrei eventualmente intercettare
### il get item per convertirli ai nomi
### ed aggiungo supporto per lettura e salvataggio da un csv tabellare
### a due dimensioni... in quel caso il dtype... vediamo che succede
### usando le funzioni numpy classiche...
###
### Questo differenzia questa classe dagli array cgh che
### usano un array numpy solo come rappresentazione interna.

### Ma se eredito che succede quando concateno o faccio altre cose strane??
### potrebbe essere una menata gestire queste eccezioni.
### Ad esempio x = x[1:] che fa? x[1:] dovrebbe ritornare una sottoclasse
### dello stesso tipo con la mappa aggiornata senza una colonna...
### e se faccio r_[x, x].... che succede?

### Di cosa ho bisogno: lettura scrittura csv, non dimenticarmi delle etichette
### larry sarebbe più che sufficiente se non fosse che legge una merda i csv.

### Forse conviene anche a me usare un array numpy e fornire un minimo di
### funzionalità nella gestione delle etichette
### Per adesso potremmo anche lasciarla perdere!!!!!!!!!!!
### Alla fine la lettura e la scrittura sono sempre dipendenti dall'applicazione
### Un po' di codice per usare un larry alla fine non è niente di complicato,
### anche dovendo leggere manualmente il csv.

### Ad esempio la classe mi è tornata comoda, ma solo per converitre un csv
### in un altro ma era poco avanzata per fare operazioni davvero furbe,
### tipo estrarre direttamente la sottomarice perché non potevo usarla
### come un ndarray... però sottoclassarlo rende le cose complicate.

## Forse ereditare è piùà facile di quel che sembra

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
        return np.copy(self._rdata)






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
