from nose.tools import *

import warnings
import numpy as np
from cghutils import ArrayCGH

class TestArrayCGH(object):

    def test_init(self):
        with warnings.catch_warnings(record=True) as w:
            acgh = ArrayCGH({'row': range(11)})
            assert issubclass(w[-1].category, RuntimeWarning)

        assert_equals(10, acgh.metadata['Row num'])
        assert_equals(0, acgh.metadata['Col num'])

        # At startup all the "columns" of the table are empty
        # it's like a clean array.
        assert_equals(0, len(acgh.ref_signal))
        assert_true(np.all(np.isnan(acgh.ref_signal)))

    def test_init_with_values(self):
        pass

    #def test_wrong_name(self):
        #assert_raises(AttributeError, getattr, self.acgh, 'foo_signal')

        #ref_signal=foo_signal,
        #sample_signal=foo_signal,
        #ref_bg=foo_signal,
        #sample_bg=foo_signal,
        #chromosomes=[1]*100,
        #start_base=foo_signal,
        #end_base=foo_signal)

        #for x in (acgh.ref_signal, acgh.sample_signal,
        #          acgh.ref_bg, acgh.sample_bg):
        #    assert_equals(100, len(x))
        #    assert_equals(foo_signal, x)
