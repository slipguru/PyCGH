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


#-----------
from .datatypes.arraycgh import ArrayCGH
from .datatypes.cytobands import _check_label, ChromosomeBand
import random as rnd
import itertools as it

def sampler(pmf):
    import bisect
    #pmf = ((0, 0.1), (1, 0.2), (2, 0.5), (3, 0.2))
    #pmf = ((2, 1.0),) # default

    pmf = sorted(pmf)
    events, prob = zip(*pmf)
    cdf = np.cumsum(prob)
    assert(cdf[-1] == 1.0)

    def _sampler(p):
        index = bisect.bisect(cdf, p)
        if p == cdf[index-1]: # bound check
            return events[index-1]
        else:
            return  events[index]

    return _sampler

    #p = np.random.random()
    #test_value = event(p)

class ArrayCGHSynth(object):

    def __init__(self, geometry, design, alterations=None, cytostructure=None):
        self._nrow, self._ncol = geometry

        design = dict(design) # ensure dict structure

        if not alterations is None and cytostructure is None:
            raise ValueError('mandatory')
        elif alterations is None:
            alterations = tuple() # default: empty iterable!

        # Fullfilled id-list (with standard "unused ids" )
        self._id = (design.keys() +
                    ['--'] * ((self._nrow * self._ncol) - len(design)))
        # Associated mask
        self._mask = np.ones(len(self._id), dtype=bool)
        self._mask[:len(design)] = False

        # Chip coordinates
        self._row, self._col = zip(*it.product(xrange(self._nrow),
                                               xrange(self._ncol)))

        # Shuffling across chip
        # sampling without replacement
        order = rnd.sample(range(len(self._id)), len(self._id))
        self._id = np.asarray([self._id[i] for i in order])
        self._mask = np.asarray([self._mask[i] for i in order])

        # Coordinates (following shuffled ids order)
        self._chr = []
        self._sb = []
        self._eb = []
        for id in self._id:
            if id in design:
                chr, sb, eb = design[id]
            else:
                chr = sb = eb = -1

            self._chr.append(chr)
            self._sb.append(sb)
            self._eb.append(eb)

        self._chr = np.asarray(self._chr)
        self._sb = np.asarray(self._sb)
        self._eb = np.asarray(self._eb)

        # Try to merge following loop with the previous one!

        # For each probe we check if it belongs to a specified
        # alteration

        samplers = dict((a, (alterations[a], [])) for a in alterations)
        for i in xrange(len(self._id)):
            c, s, e = self._chr[i], self._sb[i], self._eb[i]

            if c < 0: continue
            probe_bands = cytostructure[c][s:e]

            for a in alterations:
                chr, arm, band = _check_label(a)
                altered_band = cytostructure[chr].band(arm + band)

                if any(pb in altered_band for pb in probe_bands):
                    samplers[a][1].append(i)

        for a in alterations:
            print samplers[a]
            print self._id[samplers[a][1]]

        #print self._chr
        #print samplers


    def draw(self):
        ref_sigma = np.random.uniform(0.1, 0.2)
        reference = np.random.normal(loc=2.0, scale=ref_sigma,
                                     size=sum(~self._mask))


        import bisect
        #pmf = ((0, 0.1), (1, 0.2), (2, 0.5), (3, 0.2))

        pmf = ((2, 1.0),) # default

        pmf = sorted(pmf)
        events, prob = zip(*pmf)
        cdf = np.cumsum(prob)
        assert(cdf[-1] == 1.0)

        def event(p):
            index = bisect.bisect(cdf, p)
            if p == cdf[index-1]: # bound check
                return events[index-1]
            else:
                return  events[index]

        p = np.random.random()
        test_value = event(p)


        #test_value = 2.0
        tissue_prop = np.random.uniform(0.3, 0.7)

        test_sigma = np.random.uniform(0.1, 0.2)
        test_error = np.random.normal(loc=0.0, scale=test_sigma,
                                      size=sum(~self._mask))

        test = (np.ones_like(test_error) * (test_value * tissue_prop +
                                            2 * (1.0 - tissue_prop)) +
                test_error)

        tmp = -np.ones(self._nrow * self._ncol)

        tmpr = tmp.copy()
        tmpr[~self._mask] = reference
        reference = tmpr

        tmpt = tmp.copy()
        tmpt[~self._mask] = test
        test = tmpt


        #id, row, col, reference_signal, test_signal,
        #         chromosome, start_base, end_base, mask=None, **kwargs):
        return ArrayCGH(id = self._id,
                        row = self._row,
                        col = self._col,

                        reference_signal = reference,
                        test_signal = test,

                        chromosome = self._chr,
                        start_base = self._sb,
                        end_base = self._eb,
                        mask = self._mask)
