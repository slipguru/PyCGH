from nose.tools import *

from cghutils.utils import CytoBands

class TestCytoBands(object):

    def setup(self):
        self.cb = CytoBands()

    def test_bands(self):
        assert_equals(['p36.32', 'p36.31'], self.cb[1][2300001:7200000])
        assert_equals(['p36.33', 'p36.32', 'p36.31'], self.cb[1][2300000:7200000])
        assert_equals(['p36.33', 'p36.32', 'p36.31', 'p36.23'], self.cb[1][2300000:7200001])
        assert_equals(['q11.21'], self.cb['Y'][13500000:14000000])

    def test_bands_default(self):
        assert_equals(62, len(self.cb[1][2300001:]))
        assert_equals(3, len(self.cb[1][:7200000]))
        assert_equals(63, len(self.cb[1])) #ALL

        assert_equals(63, len(self.cb['1'])) #ALL
        assert_equals(11, len(self.cb['Y'])) #ALL

    def test_reverse_query(self):
        assert_equals((1, 249250621), self.cb[1].limits())
        assert_equals((2300001, 5400000), self.cb[1].limits('p36.32'))
        assert_equals((1, 7200000), self.cb[1].limits('p36.3'))
        assert_equals((1, 28000000), self.cb[1].limits('p36'))
        assert_equals((1, 84900000), self.cb[1].limits('p3'))
        assert_equals((1, 125000000), self.cb[1].limits('p'))

        assert_equals(self.cb[1].limits(), self.cb[1].limits(''))

        assert_raises(ValueError, self.cb[1].limits, 'f')
        assert_raises(ValueError, self.cb[1].limits, '2')
        assert_raises(ValueError, self.cb[1].limits, 2)
        assert_raises(ValueError, self.cb[1].limits, 'p7')
        assert_raises(ValueError, self.cb[1].limits, 'p3.')
        assert_raises(ValueError, self.cb[1].limits, 'p36.')

    def test_get_subbands(self):
        assert_equals(['p36.33', 'p36.32', 'p36.31'], self.cb[1].sub_bands('p36.3'))

        bands_table = self.cb[1].sub_bands('p36.3', with_limits=True)
        assert_equals(3, len(bands_table))
        assert_equals(('p36.33', 'p36.32', 'p36.31'), zip(*bands_table)[0])

        assert_raises(ValueError, self.cb[1].sub_bands, 'p36.')
        assert_raises(ValueError, self.cb[1].sub_bands, 'z36')