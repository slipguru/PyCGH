import itertools as it
import operator as op

from ..utils import _file_handle

# Cytostructure is a parser
# ChromosomeStructure is the resulting object
# ChromosomeBand is a representation of a band (for all resolution)

# Cytoband access code --------------------------------------------------------
class ChromosomeBand(object):
    def __init__(self, chromosome, label, start_base, end_base,
                 gstrand=None, sub_bands=None):

        if start_base >= end_base:
            raise ValueError('wrong band coordinates '
                             '%d-%d' % (start_base, end_base))

        self._chromosome = chromosome
        self._label = label
        self._start_base = start_base
        self._end_base = end_base

        self._gstrand = gstrand
        self._sub_bands = sub_bands

    @property
    def chromosome(self):
        return self._chromosome

    @property
    def label(self):
        return self._label

    @property
    def start_base(self):
        return self._start_base

    @property
    def end_base(self):
        return self._end_base

    @property
    def gstrand(self):
        return self._gstrand

    def expand(self):
        return self._sub_bands

    def __eq__(self, other):
        return (self.chromosome == other.chromosome and
                self.label == other.label and
                self.start_base == other.start_base and
                self.end_base == other.end_base)

    def __ne__(self, other):
        return (not self == other)

    def __cmp__(self, other):
        # Equality shortcut  (using __eq__)
        if self == other: return 0

        # Chromosome Order
        if self.chromosome != other.chromosome:
            raise RuntimeError('bands belong to different chromosomes')

        # Checking ambiguity
        if self.start_base <= other.start_base:
            first, last = self, other
        else:
            first, last = other, self

        # Overlapping bands
        if first.end_base >= last.start_base:
            raise RuntimeError('ambigous comparison between '
                               '%s%s and %s%s' % (self.chromosome, self.label,
                                                  other.chromosome, other.label))

        return cmp(self.start_base, other.start_base)

    def __contains__(self, item):
        return (item.chromosome == self.chromosome and
                self.start_base <= item.start_base <= self.end_base and
                self.start_base <= item.end_base <= self.end_base)

    def __str__(self):
        return '%s%s [%d-%d]' % (_int2chr(self.chromosome), self.label,
                                 self.start_base, self.end_base)

    def __repr__(self):
        return '%s%s' % (_int2chr(self.chromosome), self.label)


class ChromosomeStructure(object):

    def __init__(self, chromosome, starts, ends, labels, gstrands):
        if 1 <= chromosome <= 24:
            self._chromosome = chromosome
        else:
            raise ValueError('wrong chromosome number %d' % chromosome)

        # Sort bands by starting base
        sorted_bands = sorted(zip(labels, starts, ends, gstrands),
                              key=op.itemgetter(1))

        self._band_keys = dict((k[0], i) for i, k in enumerate(sorted_bands))
        self._bands = tuple(ChromosomeBand(self._chromosome, *band)
                            for band in sorted_bands)

    @property
    def chromosome(self):
        return self._chromosome

    @property
    def str_chromosome(self):
        return _int2chr(self._chromosome)

    def band(self, label=None):

        if label is None: # full chromosome
            return ChromosomeBand(self._chromosome, '', # empty label
                                  self._bands[0].start_base,
                                  self._bands[-1].end_base,
                                  sub_bands = self._bands)

        label = str(label)
        if label in self._band_keys:
            return self._bands[self._band_keys[label]] # O(1)
        else:
            if label.endswith('.'): #special case of error
                raise ValueError('wrong band representation %s' % label)

            # O(n)
            sub_bands = tuple(band for band in self._bands
                              if band.label.startswith(label))

            if sub_bands:
                return ChromosomeBand(self._chromosome, label,
                                      sub_bands[0].start_base,
                                      sub_bands[-1].end_base,
                                      sub_bands = sub_bands)
            else:
                raise ValueError('wrong band value %s' % label)

    def __getitem__(self, key):

        # To be consistent with the cytogenetic meaning, the slice
        # range in considere inclusive [start, stop] instead of [start, stop[
        # and the indexing starts from 1 instead of 0

        if not isinstance(key, slice):
            raise ValueError('requested slice object')
        if not key.step is None:
            raise ValueError('slice step value not allowed')

        min_sb = self._bands[0].start_base
        max_eb = self._bands[-1].end_base

        start_base = key.start if not key.start is None else min_sb
        end_base = key.stop if not key.stop is None else max_eb

        if start_base == 0:
            raise ValueError('base slicing starts from 1')

        #         sb                   eb
        # |    |   :  |    |       |   :    |       |
        #             *******************************   (s >= start_base)
        # **************************                    (e <= end_base)
        #             **************                    fully_contained
        #       ******                                  ( s <= start_base <= e)
        #                           *********           ( s <= end_base <= e)
        #       *****************************           RESULT

        def keep(band):
            s, e = band.start_base, band.end_base
            return (( s >= start_base and e <= end_base ) or  #fully_contained
                    ( s <= start_base <= e )              or
                    ( s <= end_base <= e ))

        return tuple(band for band in self._bands if keep(band))

    def __iter__(self):
        return iter(self._bands)

    def bands_iter(self, level=None):
        if level is None:
            return self.__iter__()

        if level <= 2:
            level += 1
        elif level > 2:
            level += 2 # dot
        else:
            raise ValueError('wrong level value')

        out = list()
        for label, sub_bands in it.groupby(self._bands,
                                           lambda b: b.label[:level]):
            sub_bands = tuple(sub_bands)
            out.append(ChromosomeBand(self._chromosome, label,
                                      sub_bands[0].start_base,
                                      sub_bands[-1].end_base,
                                      sub_bands = sub_bands))
        return iter(tuple(out))

    def __len__(self):
        return len(self._bands)

    def __str__(self):
        labels = dict()
        for k, g in it.groupby(self._bands, lambda b: b.label[0]):
            labels[k] = ' | '.join(b.label for b in g)

        return 'Chr %s < %s || %s >' % (self.str_chromosome,
                                        labels['p'], labels['q'])

class CytoStructure(object):
    """Chromosome cytogenetic structure manager.

    It is a map-like collection of ChromosomeStruture objects, each one
    containing information about a specified chromosome cytogenetic structure.
    This class act as a parser of file containing needed informations.

    File format description... TODO

    Parameters
    ----------
    cytofile : str or file
        File or filename containing chromosomes structure to read.
        If the filename extension is ``.gz`` the file is decompressed.
    format : 'ucsc' or None
        Cytogenetic file format. If 'ucsc' starting base of each cytogenetic
        band is incremented of one.

    """
    def __init__(self, cytofile, format='ucsc'):

        # file reading
        fh = _file_handle(cytofile)

        # format offset
        if format =='ucsc':
            start_offset = 1
        else:
            start_offset = 0

        # parsing bands coordinates
        bands = dict()
        for line in fh:
            chr, sb, eb, label, gs = line.strip().split()
            chr = _chr2int(chr[3:])

            if not chr in bands:
                bands[chr] = {'start': [], 'end': [],
                              'label': [], 'gstrand': []}

            bands[chr]['start'].append(int(sb) + start_offset)
            bands[chr]['end'].append(int(eb))
            bands[chr]['label'].append(label)
            bands[chr]['gstrand'].append(gs)

        # Creation of a dictionary of ChromosomeStructure objects
        self._bands = dict()
        for chr in bands:
            self._bands[chr] = ChromosomeStructure(chr,
                                               bands[chr]['start'],
                                               bands[chr]['end'],
                                               bands[chr]['label'],
                                               bands[chr]['gstrand'])
    def __getitem__(self, chr):
        try:
            chr = int(chr)
        except ValueError:
            if chr.strip() == 'X': chr = 23
            elif chr.strip() == 'Y': chr = 24

        return self._bands[chr]

    def __iter__(self):
        chromosomes = sorted(self._bands.keys())
        return iter(self._bands[chr] for chr in chromosomes)

    def __len__(self):
        return len(self._bands)

# Private utils ---------------------------------------------------------------
def _chr2int(value):
    converted = str(value).strip()
    if converted == 'X':
        return 23
    elif converted == 'Y':
        return 24
    elif converted in [str(x) for x in xrange(1, 23)]:
        return int(converted)
    else:
        raise ValueError('incorrect chromosome number %s'% converted)

def _int2chr(value):
    converted = int(value)

    if converted < 23:
        return str(value)
    elif converted == 23:
        return 'X'
    elif converted == 24:
        return 'Y'
    else:
        raise ValueError('incorrect chromosome number %d'% converted)


def _check_label(label):
    label = label.strip()

    if 'p' in label:
        chr, band = label.split('p')
        arm = 'p'
    elif 'q' in label:
        chr, band = label.split('q')
        arm = 'q'
    elif len(label) <= 2:
        chr = label.strip()
        band = arm = ''
    else:
        raise ValueError('wrong chromosome arm in %s' % label)

    if not 1 <= _chr2int(chr) <= 24:
        raise ValueError('wrong chromosome number in %s' % label)

    if band: # not empty
        if band.endswith('.'):
            raise ValueError('wrong chromosome band code in %s' % label)

        # check dot position
        if len(band) >= 3 and not band[2] == '.':
            raise ValueError('wrong chromosome band code in %s' % label)

        try:
            int(band.replace('.', ''))
        except ValueError:
            raise ValueError('wrong chromosome band code in %s' % label)

    return chr, arm, band
