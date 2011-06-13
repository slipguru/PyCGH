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
    def __init__(self, names, data_type=np.float32):
        self._names = dict((n, i) for i, n in enumerate(names))
        self._samples = {}

        self._rdata = np.empty((0, len(names)), dtype=data_type)

    def append(self, name, values):
        " ADD "
        values = np.asanyarray(values)
        values.shape = (1, len(values))

        if name in self._samples:
            raise ValueError('sample already present')

        self._rdata = np.r_[self._rdata, values]
        self._samples[name] = (len(self._rdata) - 1)

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
        " DELETE "
        if not key in self._samples:
            raise ValueError('sample not existent')

        self._rdata = np.delete(self._rdata, self._samples[key], 0)

        # Shift indici... se spostasti l'ultimo al posto dal cancellato
        # ed elimino solo una riga?
        sh_keys = [k for k in self._samples if self._samples[k] > self._samples[key]]
        for k in sh_keys:
            self._samples[k] -= 1

        del self._samples[key]

    def __len__(self):
        return len(self._rdata)
