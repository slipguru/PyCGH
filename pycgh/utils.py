import itertools as it

import numpy as np

from rpy2 import robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri

def lowess(x, y, **kwargs):
    """ Lowess from R """
    rlowess = ro.r['lowess']
    return rlowess(x, y, **kwargs)

def loess(values, x, y=None, **kwargs):
    """ Loess from R """
    loess = ro.r['loess']
    predict = ro.r['predict']

    if y is None:
        fmla = ro.Formula('v ~ x')
        env = {'v': values, 'x': x}
        df = ro.DataFrame({'x': x})
    else:
        fmla = ro.Formula('v ~ x*y')
        env = {'v': values, 'x': x, 'y': y}
        df = ro.DataFrame({'x': x, 'y': y})

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
from .datatypes.cytobands import _check_label, ChromosomeBand, _chr2int
import random as rnd
import itertools as it
from collections import defaultdict

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

    def __init__(self, geometry, design,
                 alterations=None, cytostructure=None,
                 noise=(0.1, 0.2), tissue_prop=(0.3, 0.7),
                 dyes=(100, 500), outliers=(0.2, 0.5)):

        # Check Geometry
        self._nrow, self._ncol = (int(x) for x in geometry)
        if self._nrow < 0 or self._ncol < 0:
            raise ValueError('wrong geometry (%s, %s)' % (self._nrow,
                                                          self._ncol))

        # Check Sigma
        self._Smin, self._Smax = sorted(noise)
        if self._Smin < 0 or self._Smax < 0:
            raise ValueError('wrong noise extremes (%s, %s)' % (self._Smin,
                                                                self._Smax))

        # Check Tissue Proportion
        self._Tmin, self._Tmax = sorted(tissue_prop)
        if not 0 <= self._Tmin <=1 or not 0 <= self._Tmax <= 1:
            raise ValueError('wrong tissue proportion extremes '
                             '(%s, %s)' % (self._Tmin, self._Tmax))

        # Check Dyes
        self._Dmin, self._Dmax = sorted(int(x) for x in dyes)
        if self._Dmin < 0 or self._Dmax < 0:
            raise ValueError('wrong dyes extremes (%s, %s)' % (self._Dmin,
                                                               self._Dmax))

        # Check Outliers
        self._Omin, self._Omax = sorted(outliers)
        if self._Omin < 0 or self._Omax < 0:
            raise ValueError('wrong outliers noise extremes '
                             '(%s, %s)' % (self._Omin, self._Omax))

        design = dict(design) # ensure dict structure
        CHIP_LEN = (self._nrow * self._ncol)

        if alterations and not cytostructure:
            raise ValueError('missing cytostructure reference')

        # Fullfilled id-list (with standard "unused ids" )
        self._id = np.asarray((design.keys() +
                              ['--'] * (CHIP_LEN - len(design)))) # unsused ids

        # Associated mask
        self._mask = np.ones(len(self._id), dtype=bool)
        self._mask[:len(design)] = False

        # Coordinates
        self._chr = -np.ones(len(self._id), dtype=int)
        self._sb =  -np.ones(len(self._id), dtype=int)
        self._eb =  -np.ones(len(self._id), dtype=int)

        chr, sb, eb = zip(*design.values()) # Same order as keys
        self._chr[:len(design)] = chr #[_chr2int(c) for c in chr]
        self._sb[:len(design)] = sb
        self._eb[:len(design)] = eb

        # Chip coordinates
        self._row, self._col = zip(*it.product(xrange(self._nrow),
                                               xrange(self._ncol)))

        # Shuffling order across chip (using sampling without replacement)
        order = np.asarray(rnd.sample(xrange(len(self._id)), len(self._id)))

        # For each probe we check if it belongs to a specified alteration
        # This step is perfomed befor shuffling to iterate only on
        # valid probes
        # resulting variable 'samplers' is a list of pairs
        # (sampler_fun, probe_indexes)
        if alterations:
            order_masked = order[~self._mask[order]]
            rev_order = np.argsort(order_masked)
            indexes = defaultdict(list)

            valid_clones_indexes = np.arange(len(self._id))[~self._mask]
            for i in valid_clones_indexes:
                c, s, e = self._chr[i], self._sb[i], self._eb[i]

                probe_bands = cytostructure[c][s:e]

                for a in alterations:
                    chr, arm, band = _check_label(a) # split-label
                    altered_band = cytostructure[chr].band(arm + band)

                    if any(pb in altered_band for pb in probe_bands):
                        indexes[a].append(rev_order[i])

            self._samplers = [(sampler(alterations[a]),
                               np.asarray(indexes[a], dtype=int))
                              for a in alterations]

        else:
            self._samplers = []

        # Then we can perform data shuffling across chip
        self._id = self._id[order]
        self._mask = self._mask[order]
        self._chr = self._chr[order]
        self._sb = self._sb[order]
        self._eb = self._eb[order]

    @property
    def geometry(self):
        return self._nrow, self._ncol

    @property
    def noise(self):
        return self._Smin, self._Smax

    @property
    def tissue_prop(self):
        return self._Tmin, self._Tmax

    @property
    def dyes(self):
        return self._Dmin, self._Dmax

    @property
    def outliers(self):
        return self._Omin, self._Omax

    def draw(self):
        # valid number of clones
        C = sum(~self._mask)

        # -- Base (noisy) Reference signal --
        ref_sigma = np.random.uniform(self._Smin, self._Smax)
        r = np.random.normal(2.0, ref_sigma, size=C)

        # -- Base (noisy) Test signal --
        t = 2.0 * np.ones(C)

        # * Adding alterations (resampling)
        for sampler, indexes in self._samplers:
            t[indexes] = sampler(np.random.random())

        # * Adding tissue proportion bias
        tissue_prop = np.random.uniform(self._Tmin, self._Tmax)
        test_sigma = np.random.uniform(self._Smin, self._Smax)

        t = ((t * tissue_prop) +
             (   2 * (1.0 - tissue_prop)) +
             np.random.normal(0.0, test_sigma, size=C))

        # -- Noises ---

        # * Adding outliers
        O = np.random.uniform(0.2, 0.5)
        r += np.abs(np.random.normal(0.0, O, size=C))
        t += np.abs(np.random.normal(0.0, O, size=C))

        # * Wave effect
        a = np.random.uniform(0., (np.log2(t) - np.log2(r)).std())
        kl = 2./max(self._eb[~self._mask])

        w = (a * np.sin(kl * np.pi * self._sb[~self._mask]) +
             np.random.normal(0.0, (a/2.)**2, size=C))

        r += 2**w
        t += 4**w

        # * Multiplicative Dye Bias
        r *= np.random.uniform(self._Dmin, self._Dmax)
        t *= np.random.uniform(self._Dmin, self._Dmax)

        # Thresholding!
        r[r < 0] = 0.0
        t[t < 0] = 0.0

        # -- Producing final signal --
        wave = -np.ones(len(self._id), dtype=float)
        wave[~self._mask] = w

        reference = -np.ones_like(wave)
        reference[~self._mask] = r

        test = -np.ones_like(wave)
        test[~self._mask] = t

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
                        mask = self._mask,
                        wave = wave)
