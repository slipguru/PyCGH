import numpy as np
from numpy.testing import *

from ..readers import agilent
from ..analysis import cghnormaliter
from ..datatypes import ArrayCGH

# importing the path for the small testing agilent sample
from test_agilent import SAMPLE_PATH

def _test_instantiation():
    acgh = agilent(SAMPLE_PATH)
    acgh['norm'] = cghnormaliter(acgh)

def test_duplication():
    acgh = agilent(SAMPLE_PATH)
    acgh['id'][300:350] = acgh['id'][350:400]

    idxs, norm = cghnormaliter(acgh)

    print len(acgh), acgh.size, len(norm)

    mask = np.ones(acgh.size, dtype=bool)
    mask[idxs] = False
    acgh['mask'][~acgh['mask']] = mask

    print len(acgh), acgh.size, len(norm)

    acgh['norm'] = norm

    raise Exception('TO BE CHECK!!')



