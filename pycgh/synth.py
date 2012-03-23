import random as rnd
import itertools as it
from collections import defaultdict

import numpy as np
try:
    from scipy import linalg as la
except ImportError:
    from numpy import linalg as la

from .datatypes.arraycgh import ArrayCGH
from .datatypes.cytobands import ChromosomeBand, _check_label, _chr2int

def _sampler(pmf):
    import bisect

    pmf = sorted(pmf)
    events, prob = zip(*pmf)
    cdf = np.cumsum(prob)

    def sampler(p):
        index = bisect.bisect(cdf, p)
        if p == cdf[index-1]: # bound check
            return events[index-1]
        else:
            return  events[index]

    return sampler

def _mask_signal(signal, mask,
                 missing_value=ArrayCGH.MISSING_FLOAT, dtype=float):
    full_signal = np.empty(len(mask), dtype=dtype)
    full_signal.fill(missing_value)
    full_signal[~mask] = signal

    return full_signal

def _mvnpdf(mu, cov):
    mu = np.asarray(mu)
    cov = np.asarray(cov)

    # Precalculation
    covI = la.inv(cov)

    def mvnpdf(x, y):
        dev = (np.c_[x, y] - mu)
        return np.exp(-0.5 * (np.dot(dev, covI) * dev).sum(axis=1))

    return mvnpdf

class ArrayCGHSynth(object):

    def __init__(self, geometry, design,
                 alterations=None, cytostructure=None,
                 noise=(0.1, 0.2), tissue_prop=(0.3, 0.7),
                 outliers_prop=(0.001, 0.002),
                 dyes=(100, 500), spatial_bias=0.5):

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

        # Check Outliers Proportion
        self._Omin, self._Omax = sorted(outliers_prop)
        if not 0 <= self._Omin <=1 or not 0 <= self._Omax <= 1:
            raise ValueError('wrong outliers proportion extremes '
                             '(%s, %s)' % (self._Omin, self._Omax))

        # Check Dyes
        self._Dmin, self._Dmax = sorted(int(x) for x in dyes)
        if self._Dmin < 0 or self._Dmax < 0:
            raise ValueError('wrong dyes extremes (%s, %s)' % (self._Dmin,
                                                               self._Dmax))

        # Spatial Bias Probability
        self._SB = spatial_bias
        if not 0.0 <= self._SB <= 1.0:
            raise ValueError('wrong spatial bias probability %s' % self._SB)

        design = dict(design) # ensure dict structure
        CHIP_LEN = (self._nrow * self._ncol)

        # Checking and filling alteration probabilities
        if alterations:
            if not cytostructure:
                raise ValueError('missing cytostructure reference')

            for a in alterations:
                levels, p = zip(*alterations[a])

                p_sum = sum(p)
                if p_sum > 1.0:
                    raise ValueError("sum of probabilities for "
                                     "'%s' greater than 1.0" % a)
                elif p_sum < 1.0:
                    if 2 in levels:
                        raise ValueError("sum of probabilities for "
                                         "'%s' less than 1.0" % a)
                    else: #filling
                        alterations[a].append((2, 1.0 - p_sum))

        # Fullfilled id-list (with standard "unused ids" )
        missing_num = (CHIP_LEN - len(design))
        self._id = np.asarray((design.keys() +      # unsused ids
                               [ArrayCGH.MISSING_STRING] * missing_num))

        # Associated mask
        self._mask = np.ones(len(self._id), dtype=bool)
        self._mask[:len(design)] = False

        # Coordinates
        self._chr = ArrayCGH.MISSING_INT * np.ones(len(self._id), dtype=int)
        self._sb = ArrayCGH.MISSING_INT * np.ones(len(self._id), dtype=int)
        self._eb = ArrayCGH.MISSING_INT * np.ones(len(self._id), dtype=int)

        chr, sb, eb = zip(*design.values()) # Same order as keys
        self._chr[:len(design)] = chr #[_chr2int(c) for c in chr]
        self._sb[:len(design)] = sb
        self._eb[:len(design)] = eb

        # Chip coordinates
        self._row, self._col = zip(*it.product(xrange(self._nrow),
                                               xrange(self._ncol)))
        self._row = np.asarray(self._row)
        self._col = np.asarray(self._col)

        # Shuffling order across chip (using sampling without replacement)
        order = np.asarray(rnd.sample(xrange(len(self._id)), len(self._id)))

        # For each probe we check if it belongs to a specified alteration
        # This step is perfomed iteratint only on valid probes
        # Resulting variable is a list of pairs (sampler_fun, probe_indexes)
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

            self._samplers = [(_sampler(alterations[a]),
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
    def outliers_prop(self):
        return self._Omin, self._Omax

    @property
    def dyes(self):
        return self._Dmin, self._Dmax

    @property
    def spatial_bias(self):
        return self._SB

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
        true_test_signal = _mask_signal(t, self._mask)  # Saving true Signal

        # * Adding tissue proportion bias
        tissue_prop = np.random.uniform(self._Tmin, self._Tmax)
        test_sigma = np.random.uniform(self._Smin, self._Smax)

        t = ((t * tissue_prop) +
             (   2 * (1.0 - tissue_prop)) +
             np.random.normal(0.0, test_sigma, size=C))

        # -- Noises ---

        # * Wave effect
        pos_tr = np.ma.masked_less(t/r, 0.0) # Only Positive values
        a = np.random.uniform(0., pos_tr.std()/2.)
        kl = 2./max(self._eb[~self._mask])

        w = (a * np.sin(kl * np.pi * self._sb[~self._mask]) +
             np.random.normal(0.0, (a/2.)**2, size=C))

        r += 2**w
        t += 4**w

        # * Spatial bias
        for signal in (r, t):
            if np.random.uniform(0.0, 1.0) > self._SB:
                # Trend position
                mu = np.array([np.random.randint(0, self._nrow),
                               np.random.randint(0, self._ncol)])

                # Trend shape
                var = (np.random.uniform(0, self._nrow),
                       np.random.uniform(0, self._ncol))
                cov = np.random.uniform(-1, 1) * np.mean(var)
                Sigma = np.array([[var[0], cov],
                                  [cov, var[1]]])

                # Bias calculation (random peak orientation)
                mvn = _mvnpdf(mu, Sigma)
                bias = mvn(self._row[~self._mask],
                           self._col[~self._mask]) * rnd.choice([1, -1])

                # Randomly applied on reference or test
                signal += (bias + np.random.normal(0.0, 0.1, size=C))

        # * Multiplicative Dye Bias
        r *= np.random.uniform(self._Dmin, self._Dmax)
        t *= np.random.uniform(self._Dmin, self._Dmax)

        # * Adding outliers
        for signal in (r, t):
            proportion = np.random.uniform(self._Omin, self._Omax)
            number = int(proportion * len(signal))
            indexes = rnd.sample(range(len(signal)), number)

            sigma = np.random.uniform(signal.std(), signal.std()*10.)

            signal[indexes] += np.abs(np.random.normal(0.0, sigma, size=number))

        # Thresholding!
        for signal in (r, t):
            pos = np.ma.masked_less(signal, 0.0, copy=False) # Positive values
            np.clip(signal, pos.min(), np.inf, out=signal)   # In-place

        # -- Producing final signal --
        return ArrayCGH(id = self._id,
                        row = self._row,
                        col = self._col,

                        reference_signal = _mask_signal(r, self._mask),
                        test_signal = _mask_signal(t, self._mask),

                        chromosome = self._chr,
                        start_base = self._sb,
                        end_base = self._eb,
                        mask = self._mask,

                        true_test_signal = true_test_signal)