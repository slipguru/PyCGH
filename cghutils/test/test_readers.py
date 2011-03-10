from nose.tools import *
import os

from cghutils.readers import AgilentReader

class TestAgilentReader(object):
    def setup(self):
        par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = AgilentReader(os.path.join(par_dir, 'test_agilent.txt'))

    def test_acces_params(self):
        value = self.reader.param('Protocol_Name')
        assert_equals('CGH-v4_95_Feb07 (Read Only)', value)

        value = self.reader.param('Scan_NumChannels')
        assert_equals(2, value)

        value = self.reader.param('Scan_MicronsPerPixelY')
        assert_equals(5.0, value)

        assert_raises(KeyError, self.reader.param, 'fake')

    def test_params_keys(self):
        params = self.reader.params_list()
        assert_equals(34, len(params))

    def test_access_stats(self):
        value = self.reader.stat('gDarkOffsetAverage')
        assert_equals(20.185, value)

        value = self.reader.stat('rLocalBGInlierNum')
        assert_equals(45140, value)

        assert_raises(KeyError, self.reader.stat, 'fake')

    def test_stats_keys(self):
        stats = self.reader.stats_list()
        assert_equals(143, len(stats))

    def test_access_features(self):
        value = self.reader.feature('LogRatio')
        assert_equals(4460, len(value))

        assert_raises(KeyError, self.reader.feature, 'fake')

    def test_features_keys(self):
        features = self.reader.features_list()
        assert_equals(40, len(features))


class TestGPLReader(object):
    def test_foo(self):
        pass
