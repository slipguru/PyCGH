from nose.tools import *
import os

from cghutils.readers import AgilentReader

class TestCGHReader(object):
    def setup(self):
        par_dir = os.path.split(os.path.abspath(__file__))[0]
        self.reader = AgilentReader(os.path.join(par_dir, 'test_agilent.txt'))

    def test_remove_controls(self):
        pass
