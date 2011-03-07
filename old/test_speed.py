import time
from nose import tools

import cghreader as new
import cghreader_old as old

st = time.time()
reader_new = new.CGHReader('test.txt')
print 'New', time.time() - st

st = time.time()
reader_old = old.CGHReader('test.txt')
print 'Old', time.time() - st

tools.assert_equals(reader_new.FEPARAMS,
                    reader_old.FEPARAMS)
print '- FEPARAMS OK!'

tools.assert_equals(reader_new.STATS,
                    reader_old.STATS)
print '- STATS OK!'

tools.assert_equals(sorted(reader_new.FEATURES.keys()),
                    sorted(reader_old.FEATURES.keys()))

print reader_new.FEATURES['IsManualFlag']
print reader_old.FEATURES['IsManualFlag']

for k in reader_new.FEATURES:
    tools.assert_true((reader_new.FEATURES[k] ==
                       reader_old.FEATURES[k]).all())
    print "- FEATURES['%s']\t Ok!" % k
