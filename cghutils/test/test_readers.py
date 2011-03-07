from nose.tools import *
import os

from cghutils.readers import AgilentReader

class TestCGHReader(object):
    def setup(self):
        par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = AgilentReader(os.path.join(par_dir, 'test_agilent.txt'))

    def test_acces_FEPARAMS(self):
        value = self.reader.param('Protocol_Name')
        assert_equals('CGH-v4_95_Feb07 (Read Only)', value)

        value = self.reader.param('Scan_NumChannels')
        assert_equals(2, value)

        value = self.reader.param('Scan_MicronsPerPixelY')
        assert_equals(5.0, value)

        assert_raises(KeyError, self.reader.param, 'fake')

        all_params = self.reader.param()
        assert_equals(34, len(all_params))
        assert_equals(self.reader.param('Protocol_Name'),
                      all_params['Protocol_Name'])

    def test_access_STATS(self):
        value = self.reader.stat('gDarkOffsetAverage')
        assert_equals(20.185, value)

        value = self.reader.stat('rLocalBGInlierNum')
        assert_equals(45140, value)

        assert_raises(KeyError, self.reader.stat, 'fake')

        all_stats = self.reader.stat()
        assert_equals(143, len(all_stats))
        assert_equals(self.reader.stat('gDarkOffsetAverage'),
                      all_stats['gDarkOffsetAverage'])

    def test_access_FEATURES(self):
        value = self.reader.feature('LogRatio')
        assert_equals(4460, len(value))

        #assert_raises(KeyError, self.reader.FEATURES.__getitem__, 'fake')
        #assert_equals(40, len(self.reader.FEATURES))
