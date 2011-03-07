from nose.tools import *

from agilentreader import AgilentReader

class TestCGHReader(object):
    def setup(self):
        self.reader = AgilentReader('test_agilent.txt')

    def test_acces_FEPARAMS(self):
        value = self.reader.FEPARAMS['Protocol_Name']
        assert_equals('CGH-v4_95_Feb07 (Read Only)', value)

        value = self.reader.FEPARAMS['Scan_NumChannels']
        assert_equals(2, value)

        value = self.reader.FEPARAMS['Scan_MicronsPerPixelY']
        assert_equals(5.0, value)

        assert_raises(KeyError, self.reader.STATS.__getitem__, 'fake')

        assert_equals(34, len(self.reader.FEPARAMS))

    def test_access_STATS(self):
        value = self.reader.STATS['gDarkOffsetAverage']
        assert_equals(20.185, value)

        value = self.reader.STATS['rLocalBGInlierNum']
        assert_equals(45140, value)

        assert_raises(KeyError, self.reader.STATS.__getitem__, 'fake')

        assert_equals(143, len(self.reader.STATS))

    def test_access_FEATURES(self):
        value = self.reader.FEATURES['LogRatio']
        assert_equals(4460, len(value))

        assert_raises(KeyError, self.reader.FEATURES.__getitem__, 'fake')

        assert_equals(40, len(self.reader.FEATURES))
