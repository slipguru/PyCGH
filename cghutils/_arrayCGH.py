import numpy as np
from numpy import ma

# This class represent only a useful way to read and manipulate
# the values associated with an aCGH

# The functionalities of the class are limited to reading and filtering data
# Is possible to serialize/deserialize the class in a readble tab separated
# format

# Moreover, giving a reader class is possible to initialize
# the class from other formats (e.g. Agilent Scanner output file)

# The class automatically deal with missing data using numpy masked array

# The standard constructor accept a dictionary with all the needed values
# The class was thought as a read-only structure

TYPES = [('row', np.int16),             # no missing
         ('col', np.int16),             # //
         ('id', np.str_),
         ('ref_signal', np.float_),
         ('sample_signal', np.float_),
         ('ref_bg', np.float_),
         ('sample_bg', np.float_),
         ('chromosome', np.int16),
         ('start_base', np.int_),
         ('end_base', np.int_)]

class ArrayCGH(object):

    def __init__(self, input_dict, filter=None):
        # Ordered list to record array
        names = [k[0] for k in TYPES]
        buff = [np.asarray(input_dict[k[0]], dtype=k[1]) for k in TYPES]

        try:
            self._rdata = np.rec.fromarrays(buff, names=names)
        except ValueError, e:
            msg, index = str(e).split(' array', 1)
            raise ValueError("%s '%s'" % (msg, names[int(index)]))

        if filter is None:
            filter = np.ones(len(self._rdata['row']), dtype=bool)
        elif len(filter) != len(self._rdata['row']):
            raise ValueError("invalid mask length")
        self._filter = filter

        # Utile..... se abbiamo solo una unica maschera
        # che ce ne facciamo...
        # Ma mi interessa davvero perdere tutto sto tempo?!?!
        # Se i dati li passa l'utente li passa completi:
        # controllo se row e col ci sono tutti altrimenti riempio
        # e creo la mashera (come add in R).
        # ritorno sempre delle copie per non modificare l'array
        # In lettura da file faccio lo stesso, il filling e' solo in presenza
        # di righe e colonne mancanti, l'utenti non se ne frega dei missing
        # quelli ce li metto io quando e se servono.
        # i filtri si fanno con dei metodi che al massimo espendono
        # la maschera di filtro... basta cosi'
        # Per R ho dei metodi che mi estraggono i vettori completi
        # dove l'utente puo' decidere come fare il fill
        # [key] -> ritorna filtrati
        # filled(key, fill_values=XX, dtype=Casting) -> fillato....
        self._rdata = np.ma.masked_array(np.asarray(self._rdata), mask=~filter)

        # possibles optimizations:
        #   self._rdata = np.asarray(self._rdata)
        #   self._rdata = self._rdata.view(ndarray)

    # the metod returns a copy
    def __getitem__(self, key):
        return self._rdata[key]#[self._filter]

    def compute_ratios(self, bg_subtracted=False):
        ref_signal = self._rdata['ref_signal']
        sample_signal = self._rdata['sample_signal']

        if bg_subtracted:
            ref_signal -= self._rdata['ref_bg']
            sample_signal -= self._rdata['sample_bg']

        return ma.log2(sample_signal / ref_signal)

    #def __init__(self, input_dict):
    #
    #    # Set of expected values
    #    acgh_keys = set(('row', 'col',
    #                     'ref_signal', 'sample_signal',
    #                     'ref_bg', 'sample_bg',
    #                     'chromosome', 'start_base', 'end_base'))
    #
    #    # Check if there exist at least one valid list of values
    #    valid_keys = list(acgh_keys.intersection(input_dict))
    #    if not valid_keys:
    #        raise ValueError('input dictionary does not contain expected keys')
    #
    #    # Check if some lists of values are missing (others will be discarded)
    #    missing = acgh_keys.difference(input_dict)
    #    if missing:
    #        warnings.warn('missing lists of values: %s' % ', '.join(missing),
    #                      RuntimeWarning)
    #
    #    # Array Geometry
    #    row_num = np.max(input_dict.get('row', [0]))
    #    col_num = np.max(input_dict.get('col', [0]))
    #    array_length = (row_num*col_num) or len(input_dict[valid_keys[0]])
    #
    #    # Data population
    #    empty_array = np.array([np.nan]*(array_length))
    #    self.metadata = {'Row num': row_num,
    #                     'Col num': col_num}
    #
    #    for key in valid_keys:
    #        setattr(self, key, np.asarray(input_dict[key], dtype=float))
    #
    #    for key in missing:
    #        setattr(self, key, self._empty_array)
    #
    #    # Masked arrays
    #    #a = np.array([1, 2, np.nan])
    #    #a = [1, 2, None]
    #    #b = ma.masked_invalid(np.asarray(a, dtype=float))
    #
