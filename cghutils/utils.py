import itertools as it

import numpy as np

from rpy2 import robjects
import rpy2.robjects.numpy2ri

def lowess(x, y, **kwargs):
    """ Lowess from R """
    rlowess = robjects.r['lowess']
    return rlowess(x, y, **kwargs)

def loess(values, x, y=None, **kwargs):
    """ Loess from R """
    loess = robjects.r['loess']
    predict = robjects.r['predict']

    if y is None:
        fmla = robjects.Formula('v ~ x')
        env = {'v': values, 'x': x}
        df = robjects.DataFrame({'x': x})
    else:
        fmla = robjects.Formula('v ~ x*y')
        env = {'v': values, 'x': x, 'y': y}
        df = robjects.DataFrame({'x': x, 'y': y})

    for k in env:
        fmla.environment[k] = env[k]

    out = loess(fmla, **kwargs)
    return np.asarray(predict(out, df))

def array_trend(values, col, row):
    return loess(values, col, row, span=0.03, degree=1,
                 normalize=True, family='gaussian', iterations=3)

def probes_average(probes_id, probes_values, avg_function=np.mean):
    summary = dict()
    for id, v in it.izip(probes_id, probes_values):
        summary.setdefault(id, []).append(v)

    return dict((id, avg_function(summary[id])) for id in summary)


# Cytoband access code --------------------------------------------------------
class ChromosomeBands(object):
    def __init__(self, chromosome, starts, ends, names):
        self._chromosome = chromosome
        self._starts = starts
        self._ends = ends
        self._names = names

    def __getitem__(self, key):

        # To be consistent with the cytogenetic meaning, the slice
        # range in considere inclusive [start, stop] instead of [start, stop[
        # and the indexing starts from 1 instead of 0

        if not isinstance(key, slice):
            raise ValueError('requested slice object')
        if not key.step is None:
            raise ValueError('slice step value not allowed')

        start_base = key.start if not key.start is None else min(self._starts)
        end_base = key.stop if not key.stop is None else max(self._ends)

        #         sb                   eb
        # |    |   :  |    |       |   :    |       |
        #             *******************************   (s >= start_base)
        # **************************                    (e <= end_base)
        #             **************                    fully_contained
        #       ******                                  ( s <= start_base <= e)
        #                           *********           ( s <= end_base <= e)
        #       *****************************           RESULT

        map = (( s >= start_base and e <= end_base ) or  #fully_contained
               ( s <= start_base <= e )              or
               ( s <= end_base <= e )
                        for s,e in zip(self._starts, self._ends))
        return [name for name, m in zip(self._names, map) if m]

    def limits(self, band=None):
        if not band:
            return min(self._starts), max(self._ends)

        try:
            band = str(band)
            if band.endswith('.'): raise ValueError() #special case of error

            min_limit = min(s for s, n in it.izip(self._starts, self._names)
                                          if n.startswith(str(band)))
            max_limit = max(e for e, n in it.izip(self._ends, self._names)
                                          if n.startswith(str(band)))
        except ValueError:
            raise ValueError('wrong band representation')

        return min_limit, max_limit

    def sub_bands(self, band='', with_limits=False):
        band = str(band)
        if band.endswith('.'):
            raise ValueError('wrong band representation') #special case of error

        if with_limits:
            table = it.izip(self._names, self._starts, self._ends)
            out = [(n, s, e) for n, s, e in table if n.startswith(str(band))]
        else:
            out = [n for n in self._names if n.startswith(str(band))]

        if not out:
            raise ValueError('wrong band representation') # no bands found
        return out

    def __len__(self):
        return len(self._names)

    def __str__(self):
        return str(self._names)


class CytoBands(object):
    @staticmethod
    def _chr2int(value):
            converted = value[3:]
            if converted == 'X': return 23
            if converted == 'Y': return 24
            return int(converted)

    def __init__(self, release='hg19'):
        from cghutils import UCSC

        if not release in UCSC:
            raise ValueError('wrong release specified')
        file_path = UCSC[release]['cytoBand']

        # Parsing and bands reading
        bands = dict()
        for line in open(file_path):
            chr, sb, eb, name, gs = line.strip().split('\t')
            chr = CytoBands._chr2int(chr)

            if not chr in bands:
                bands[chr] = {'start': [], 'end': [],  'name': []}

            bands[chr]['start'].append(int(sb)+1)
            bands[chr]['end'].append(int(eb))
            bands[chr]['name'].append(name)

        # Dictionary of ChromosomeBands objects creation
        self._bands = dict()
        for chr in bands:
            self._bands[chr] = ChromosomeBands(chr,
                                               bands[chr]['start'],
                                               bands[chr]['end'],
                                               bands[chr]['name'])
    def __getitem__(self, key):
        try:
            key = int(key)
        except ValueError:
            if key.strip() == 'X': key = 23
            elif key.strip() == 'Y': key = 24

        return self._bands[key]


#------------------------------------------------------------------------------
class LabeledMatrix(object):
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
