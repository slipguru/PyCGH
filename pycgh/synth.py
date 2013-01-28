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

def _sorted_pair(value):
    #operator.isNumberType(obj) py2.6
    import numbers
    if isinstance(value, numbers.Number):
        return (value, value)
    return sorted(value)

class ArrayCGHSynth(object):

    def __init__(self, geometry, design,
                 alterations=None,
                 cytostructure=None,
                 tissue_proportion=(0.3, 0.7),
                 spatial_bias_probability=0.5,
                 wave_bias_amplitude=(0.0, 0.025),
                 dye_intensity=(50, 300),
                 noise=(0.15, 0.35),
                 random_state=None):

        # Check Geometry
        self._nrow, self._ncol = (int(x) for x in geometry)
        if self._nrow < 0 or self._ncol < 0:
            raise ValueError('wrong geometry (%s, %s)' % (self._nrow,
                                                          self._ncol))

        # Check Tissue Proportion
        self._Tmin, self._Tmax = _sorted_pair(tissue_proportion)
        if not 0 <= self._Tmin <= 1 or not 0 <= self._Tmax <= 1:
            raise ValueError('wrong tissue proportion extremes '
                             '(%s, %s)' % (self._Tmin, self._Tmax))

        # Check Spatial Bias Probability
        self._SBP = spatial_bias_probability
        if not 0.0 <= self._SBP <= 1.0:
            raise ValueError('wrong spatial bias probability %s' % self._SBP)

        # Check Wave Bias Amplitude
        self._Wmin, self._Wmax = _sorted_pair(wave_bias_amplitude)
        if self._Wmin < 0 or self._Wmax < 0:
            raise ValueError('wrong wave amplitude extremes '
                             '(%s, %s)' % (self._Wmin, self._Wmax))

        # Check Dyes
        self._Dmin, self._Dmax = (int(x) for x in _sorted_pair(dye_intensity))
        if self._Dmin < 0 or self._Dmax < 0:
            raise ValueError('wrong dyes extremes (%s, %s)' % (self._Dmin,
                                                               self._Dmax))

        # Check Noise
        self._Nmin, self._Nmax = _sorted_pair(noise)
        if self._Nmin < 0 or self._Nmax < 0:
            raise ValueError('wrong noise extremes (%s, %s)' % (self._Nmin,
                                                                self._Nmax))

        # Creating Random State
        if random_state is None:
            self._rnd = np.random.mtrand._rand # already global random state
        elif isinstance(random_state, np.random.RandomState):
            self._rnd = random_state
        else:
            self._rnd = np.random.RandomState(random_state)

        design = dict(design) # ensure dict structure
        CHIP_LEN = (self._nrow * self._ncol)

        # Checking and filling alteration probabilities
        if alterations:
            for k in alterations.keys():
                if k.startswith('X') or k.startswith('Y'):
                    raise ValueError('alterations on allosomes are not '
                                     'yet allowed')

            if not cytostructure:
                raise ValueError('missing cytostructure reference')

            for a in alterations:
                try:
                    levels, p = zip(*alterations[a])
                except TypeError:
                    try:
                        alterations[a] = [(int(alterations[a]), 1.0)]
                        levels, p = zip(*alterations[a])
                    except TypeError:
                        raise ValueError('not valid list of '
                                         'alterations probability '
                                         'pairs: %s' % str(alterations[a]))

                p_sum = sum(p)
                if p_sum > 1.0:
                    raise ValueError("sum of probabilities for "
                                     "'%s' greater than 1.0" % a)
                elif p_sum < 1.0:
                    # DIFFERENCE on X and Y by gender
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
        order = np.arange(len(self._id))
        self._rnd.shuffle(order)

        # For each probe we check if it belongs to a specified alteration
        # This step is perfomed iterating only on valid probes
        # Resulting variable is a list of pairs (sampler_fun, probe_indexes)
        if alterations:
            order_masked = order[~self._mask[order]]
            rev_order = np.argsort(order_masked)
            indexes = defaultdict(list)

            valid_clones_indexes = np.arange(len(self._id))[~self._mask]
            for i in valid_clones_indexes:
                c, s, e = self._chr[i], self._sb[i], self._eb[i]

                try:
                    probe_bands = cytostructure[c][s:e]
                except KeyError:
                    raise ValueError('cytostructure cannot include all clones '
                                     'chromosome location')

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

        # Assuming Male
        self._Xindexes = (self._chr == 23)[~self._mask]
        self._Yindexes = (self._chr == 24)[~self._mask]

    @property
    def geometry(self):
        return self._nrow, self._ncol

    @property
    def noise(self):
        return self._Nmin, self._Nmax

    @property
    def tissue_proportion(self):
        return self._Tmin, self._Tmax

    @property
    def dye_intensity(self):
        return self._Dmin, self._Dmax

    @property
    def spatial_bias_probability(self):
        return self._SBP

    @property
    def wave_bias_amplitude(self):
        return self._Wmin, self._Wmax

    def draw(self, gender='male'):
        # valid number of clones
        C = sum(~self._mask)

        # -- Base signals --
        gender = gender.lower()
        r = 2.0 * np.ones(C)
        if gender == 'male':
            r[self._Xindexes] = 1.0
            r[self._Yindexes] = 1.0
        elif gender == 'female':
            r[self._Yindexes] = 0.0
        else:
            raise ValueError('not recognized gender %s' % self._gender)
        t = r.copy()

        # * Adding alterations (resampling)
        for sampler, indexes in self._samplers:
            t[indexes] = sampler(self._rnd.rand())

        # Saving true Signals
        true_test_signal = _mask_signal(t, self._mask)
        true_reference_signal = _mask_signal(r, self._mask)

        # * Adding tissue proportion bias
        tissue_prop = self._rnd.uniform(self._Tmin, self._Tmax)
        t = ((t * tissue_prop) + (r * (1.0 - tissue_prop)))

        # * Spatial bias
        for signal in (r, t):
            if self._rnd.uniform(0.0, 1.0) < self._SBP:
                # Trend position
                mu = np.array([self._rnd.randint(0, self._nrow),
                               self._rnd.randint(0, self._ncol)])

                # Shape
                vars = (self._rnd.uniform(0, (self._nrow * 2.)),
                        self._rnd.uniform(0, (self._ncol * 2.)))

                # Rotation
                theta = self._rnd.uniform(0.0, 2 * np.pi)
                U = np.array([[np.cos(theta), np.sin(theta)],
                              [-np.sin(theta), np.cos(theta)]])

                # Bias calculation (random peak orientation)
                Sigma = np.dot(np.dot(U, np.diag(vars)), U.T)
                mvn = _mvnpdf(mu, Sigma)
                bias = mvn(self._row[~self._mask],
                           self._col[~self._mask])
                scaling_factor = signal / signal.max()

                # Random intensity (bias proportional to signal intensity)
                signal += (self._rnd.uniform(-1.0, 1.0) * bias * scaling_factor)

        # * Wave effect
        a = self._rnd.uniform(self._Wmin, self._Wmax)
        freq = 8./max(self._eb[~self._mask])
        w = a * np.sin(freq * np.pi * self._sb[~self._mask])

        r *= 2**w
        t *= 4**w

        # * Signal intensity (Dye Bias + Noise)
        r_dye = self._rnd.uniform(self._Dmin, self._Dmax)
        t_dye = np.abs(r_dye + self._rnd.uniform(-r_dye/3., r_dye/3.))
        sigma = self._rnd.uniform(self._Nmin, self._Nmax)

        # Hybridization noise
        noise = sigma * np.sqrt(2.) * 0.5   # COMMENT
        response_bias = self._rnd.normal(0.0, 1.0, size=C) # Fixed?

        r *= (2.**(self._rnd.normal(np.log2(r_dye), noise, size=C) +
                   response_bias))  # Systematic Hybridization
        t *= (2.**(self._rnd.normal(np.log2(t_dye), noise, size=C) +
                   response_bias))  # Systematic Hybridization

        # * Masking not valid probes (out of signal range)
        #   All Y probes will be marked as not valid if female sample
        reference_signal = _mask_signal(r, self._mask)
        test_signal = _mask_signal(t, self._mask)
        if gender == 'female':
            invalid = _mask_signal(self._Yindexes, self._mask,
                                   False, dtype=bool)
            self._mask[invalid] = True
        self._mask[reference_signal < 0.0] = True
        self._mask[test_signal < 0.0] = True

        # -- Producing final signal --
        return ArrayCGH(id = self._id,
                        row = self._row,
                        col = self._col,

                        reference_signal = reference_signal,
                        test_signal = test_signal,

                        chromosome = self._chr,
                        start_base = self._sb,
                        end_base = self._eb,
                        mask = self._mask,

                        true_test_signal = true_test_signal,
                        true_reference_signal = true_reference_signal)
