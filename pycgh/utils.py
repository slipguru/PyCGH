import itertools as it

import numpy as np

from rpy2 import robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri

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

    def __init__(self, geometry, design, alterations=None, cytostructure=None):
        self._nrow, self._ncol = geometry

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
        order = rnd.sample(xrange(len(self._id)), len(self._id))


        # For each probe we check if it belongs to a specified alteration
        # This step is perfomed befor shuffling to iterate only on
        # valid probes
        # resulting variable 'samplers' is a list of pairs
        # (sampler_fun, probe_indexes)
        if alterations:
            rev_order = np.argsort(np.array(order))[:len(design)]
            indexes = defaultdict(list)

            for i in xrange(len(design)):
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


    def draw(self):
        # -- Reference signal --
        Smin, Smax = 0.1, 0.2 #########
        ref_sigma = np.random.uniform(Smin, Smax) #########
        reference = -np.ones(len(self._id))
        reference[~self._mask] = np.random.normal(loc=2.0, scale=ref_sigma,
                                                  size=sum(~self._mask))

        # -- Test signal --
        test = -np.ones(len(self._id))
        test[~self._mask] = 2.0 # Standard value

        # * Adding alterations (resampling)
        for sampler, indexes in self._samplers:
            test[indexes] = sampler(np.random.random())

        #print
        #print reference[~self._mask]
        #print test[~self._mask]

        # * Adding tissue proportion bias
        Tmin, Tmax = 0.3, 0.7 #########
        tissue_prop = np.random.uniform(Tmin, Tmax) #########
        test_sigma = np.random.uniform(Smin, Smax) #########

        test[~self._mask] = ((test[~self._mask] * tissue_prop) +
                             (               2  * (1.0 - tissue_prop)) +
                             np.random.normal(loc=0.0,              # Error
                                              scale=test_sigma,
                                              size=sum(~self._mask)),
                            )

        #print
        #print reference[~self._mask]
        #print test[~self._mask]

        # * Adding outliers and Dye Bias
        Bmin, Bmax = 100, 500 #########
        #reference[~self._mask] = ((reference[~self._mask] +
        #                           np.abs(np.random.normal(loc=0.0,
        #                                            scale=4, #########
        #                                            size=sum(~self._mask)))) *
        #                           np.random.uniform(Bmin, Bmax) #########
        #                         )
        #
        #test[~self._mask] = ((test[~self._mask] +
        #                      np.abs(np.random.normal(loc=0.0,
        #                                       scale=4, #########
        #                                       size=sum(~self._mask)))) *
        #                      np.random.uniform(Bmin, Bmax) #########
        #                    )

        # * Wave effect
        a = np.random.uniform(0., .5)
        kl = 3./max(self._eb) ############
        noise = np.random.normal(0.0, (a/10)**2, size=sum(~self._mask))

        w = a * np.sin(kl * np.pi * self._sb[~self._mask])

        O = np.random.uniform(0.1, 2.0)

        reference[~self._mask] = ((reference[~self._mask] +
                                   np.abs(np.random.normal(loc=0.0,
                                                    scale=O, #########
                                                    size=sum(~self._mask))) +
                                   2**(w + noise)
                                   ) *
                                   np.random.uniform(Bmin, Bmax) #########
                                 )

        test[~self._mask] = ((test[~self._mask] +
                              np.abs(np.random.normal(loc=0.0,
                                               scale=O, #########
                                               size=sum(~self._mask))) +
                              4**(w + noise)
                              ) *
                              np.random.uniform(Bmin, Bmax) #########
                            )

    wave = -np.ones_like(reference)
    wave[~self._mask] = w

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
