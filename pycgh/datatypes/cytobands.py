import itertools as it
import operator as op

# Implementation is simple but not really efficient.
# Iterations run in linear time, but to get information
# of a specific band it is a second linear time operation.

# A little improvement may be the return of ChromosomeBand objects
# instead of strings for query and iteration, to get directly
# positions and gstrand informations (if available)
class ChromosomeBand(object):
    def __init__(self, label, start_base, end_base, gstrand=None):
        self.label = label
        self.start_base = start_base
        self.end_base = end_base
        self.gstrand = gstrand
        
    # we can implement ordering methods and "startswith" methods
    # in order to reuse already implemented code
    # but simplifying user-interaction!!
    # Keeping linear time for iteration but removing position and gstrand
    # queryies
    
    # From 14:30 to 17:15

# Cytoband access code --------------------------------------------------------
class ChromosomeStructure(object):

    def __init__(self, chromosome, starts, ends, labels, gstrands):
        if 1 <= chromosome <= 24:
            self._chromosome = chromosome
        else:
            raise ValueError('wrong chromosome number %d' % chromosome)
        
        # Sort bands by starting base
        sorted_bands = sorted(zip(labels, starts, ends, gstrands),
                              key=op.itemgetter(1))
        self._labels, self._starts, self._ends, self._gstrands = zip(*sorted_bands)        
        
    @property
    def chromosome(self):
        return self._chromosome
    
    @property
    def str_chromosome(self):
        if self._chromosome < 23:
            return str(self._chromosome)
        elif self._chromosome == 23:
            return 'X'
        elif self._chromosome == 24:
            return 'Y'
        
    def __getitem__(self, key):
        
        # O(n) operation! Could be better implemented,
        # but probably this class is not a bootleneck, having
        # at maximum less than 1000 items.

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
    
        return [label for label, m in zip(self._labels, map) if m]

    def position(self, band=None):
        """ Given a band (at any resolution), returns its position on chromosome
        in terms of starting and ending bases"""
        if not band:
            return min(self._starts), max(self._ends)

        try:
            band = str(band)
            if band.endswith('.'): raise ValueError() #special case of error

            min_limit = min(s for s, l in it.izip(self._starts, self._labels)
                                          if l.startswith(str(band)))
            max_limit = max(e for e, l in it.izip(self._ends, self._labels)
                                          if l.startswith(str(band)))
        except ValueError:
            raise ValueError('wrong band representation')

        return min_limit, max_limit
    
    def gstrand(self, band):
        return self._gstrands[self._labels.index(band)]

    def expand(self, band, with_position=False):
        band = str(band)
        if band.endswith('.'):
            raise ValueError('wrong band representation') #special case of error

        if with_position:
            table = it.izip(self._labels, self._starts, self._ends)
            out = [(l, s, e) for l, s, e in table if l.startswith(str(band))]
        else:
            out = [l for l in self._labels if l.startswith(str(band))]

        if not out:
            raise ValueError('wrong band representation') # no bands found
        return out
    
    def __iter__(self):
        return iter(self._labels)
        
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
        for l in self._labels:
            label = l[:level]
            if not label in out:
                out.append(label)
        
        return iter(out)

    def __len__(self):
        return len(self._labels)

    def __str__(self):        
        plabels = [l for l in self._labels if l.startswith('p')]
        qlabels = [l for l in self._labels if l.startswith('q')]
        return 'Chr %s < %s || %s >' % (self.str_chromosome,
                                        ' | '.join(plabels),
                                        ' | '.join(qlabels),)


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
        try:
            if _is_string_like(cytofile):
                if cytofile.endswith('.gz'):
                    import gzip
                    fh = gzip.open(cytofile)
                else:
                    fh = open(cytofile, 'U')
            else:
                fh = cytofile
        except TypeError:
            raise ValueError('cytofile must be a string or file handle')

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
    converted = value.strip()
    if converted == 'X': return 23
    if converted == 'Y': return 24
    return int(converted)
    
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
        band = arm = None
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
    
    return label

def _is_string_like(obj):
    """
    Check whether obj behaves like a string.
    """
    try:
        obj + ''
    except (TypeError, ValueError):
        return False
    return True