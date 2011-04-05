import warnings
import numpy as np

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


class ArrayCGH(object):
    def __init__(self, input_dict):

        # Set of expected values
        acgh_keys = set(('row', 'col',
                         'ref_signal', 'sample_signal',
                         'ref_bg', 'sample_bg',
                         'chromosome', 'start_base', 'end_base'))

        # Check if there exist at least one valid list of values
        valid_keys = list(acgh_keys.intersection(input_dict))
        if not valid_keys:
            raise ValueError('input dictionary does not contain expected keys')

        # Check if some lists of values are missing (others will be discarded)
        missing = acgh_keys.difference(input_dict)
        if missing:
            warnings.warn('missing lists of values: %s' % ', '.join(missing),
                          RuntimeWarning)

        # Array Geometry
        row_num = np.max(input_dict.get('row', [0]))
        col_num = np.max(input_dict.get('col', [0]))
        array_length = (row_num*col_num) or len(input_dict[valid_keys[0]])

        # Data population
        empty_array = np.array([np.nan]*(array_length))
        self.metadata = {'Row num': row_num,
                         'Col num': col_num}

        for key in valid_keys:
            setattr(self, key, np.asarray(input_dict[key], dtype=float))

        for key in missing:
            setattr(self, key, self._empty_array)

        # Masked arrays
        #a = np.array([1, 2, np.nan])
        #a = [1, 2, None]
        #b = ma.masked_invalid(np.asarray(a, dtype=float))

    def __getattr__(self, key):
        try:
            return []#self._internals[key]
        except KeyError, e:
            raise AttributeError(e)
