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
        from pycgh import UCSC

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