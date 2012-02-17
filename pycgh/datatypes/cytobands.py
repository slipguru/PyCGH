import itertools as it
import operator as op

# Cytoband access code --------------------------------------------------------
class ChromosomeStructure(object):

    def __init__(self, chromosome, starts, ends, names):        
        self._chromosome = chromosome
        
        # Sort bands by starting base
        sorted_bands = sorted(zip(names, starts, ends), key=op.itemgetter(1))
        self._names, self._starts, self._ends = zip(*sorted_bands)
        
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

    def position(self, band=None):
        """ Given a band (at any resolution), returns its position on chromosome
        in terms of starting and ending bases"""
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

    def expand(self, band, with_position=False):
        band = str(band)
        if band.endswith('.'):
            raise ValueError('wrong band representation') #special case of error

        if with_position:
            table = it.izip(self._names, self._starts, self._ends)
            out = [(n, s, e) for n, s, e in table if n.startswith(str(band))]
        else:
            out = [n for n in self._names if n.startswith(str(band))]

        if not out:
            raise ValueError('wrong band representation') # no bands found
        return out
    
    def __iter__(self):
        return iter(self._names)
        
    def bands_iter(self, level=None):
        if level is None:
            return self.__iter__()
        
        if level <= 2:
            level += 1
        elif level > 2:
            level += 2
        else:
            raise ValueError('wrong level value')
        
        out = list()
        for n in self._names:
            name = n[:level]
            if not name in out:
                out.append(name)
        
        return iter(out)

    def __len__(self):
        return len(self._names)

    def __str__(self):
        return str(self._names)


class CytoStructure(object):
    
    def __init__(self, cytofile, format='ucsc'):        
        # Reading File        
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
            raise ValueError('fname must be a string, file handle, or generator')

        if format =='ucsc':
            start_offset = 1
        else:
            start_offset = 0

        # Parsing and bands reading
        bands = dict()
        for line in fh:
            chr, sb, eb, label, gs = line.strip().split()
            chr = self._chr2int(chr)

            if not chr in bands:
                bands[chr] = {'start': [], 'end': [],  'label': []}

            bands[chr]['start'].append(int(sb) + start_offset)
            bands[chr]['end'].append(int(eb))
            bands[chr]['label'].append(label)

        # Dictionary of ChromosomeBands objects creation
        self._bands = dict()
        for chr in bands:
            self._bands[chr] = ChromosomeStructure(chr,
                                               bands[chr]['start'],
                                               bands[chr]['end'],
                                               bands[chr]['label'])
    def __getitem__(self, chr):
        try:
            chr = int(chr)
        except ValueError:
            if chr.strip() == 'X': chr = 23
            elif chr.strip() == 'Y': chr = 24

        return self._bands[chr]
        
    def __len__(self):
        return len(self._bands)
        
    @staticmethod
    def _chr2int(value):
        converted = value[3:]
        if converted == 'X': return 23
        if converted == 'Y': return 24
        return int(converted)

def _is_string_like(obj):
    """
    Check whether obj behaves like a string.
    """
    try:
        obj + ''
    except (TypeError, ValueError):
        return False
    return True