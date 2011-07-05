#------------------------------------------------------------------------------
class Dataset(object):
    def __init__(self, names, data=None, sample_names=None, data_type=np.float32):
        self._names = dict((n, i) for i, n in enumerate(names))

        if data is None:
            self._rdata = np.empty((0, len(names)), dtype=data_type)
            self._samples = {}
        else:
            data = np.asanyarray(data, dtype=data_type)
            r, c = data.shape
            if len(names) != c:
                raise ValueError('data column number must be equal to number of column names')

            self._rdata = data
            if sample_names is None:
                self._samples = dict(('sample%d' % i, i) for i in xrange(1, r+1))
            else:
                self._samples = dict((name, i) for i, name in enumerate(sample_names))

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
                             sample_names=samples)

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