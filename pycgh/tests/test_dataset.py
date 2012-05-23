# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np
from numpy.testing.utils import *
from cStringIO import StringIO


from ..datatypes import DataTable
from .. import PyCGHException

DATA_TEXT = """\
# This is a Comment
*, c0, c1, c2
r0, 1, 2, 3
#########################
# This is another Comment
#########################
r1, 4, 5, 6

"""

class TestDataTable(object):
    def setup(self):
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.rlabels = ['r0', 'r1']
        self.clabels = ['c0', 'c1', 'c2']

    def test_creation(self):
        dt = DataTable(self.data, self.rlabels, self.clabels)
        assert_equal(2, dt.rnum)
        assert_equal(3, dt.cnum)
        assert_equal(2, len(dt))

        assert_equal(self.rlabels, dt.rlabels)
        assert_equal(self.clabels, dt.clabels)

    def test_labels(self):
        dt = DataTable(self.data)
        assert_equal(self.rlabels, dt.rlabels)
        assert_equal(self.clabels, dt.clabels)

    def test_labels_to_str(self):
        dt = DataTable(self.data, clabels=['col1', u'col2', 3])

        assert_equal(dt.clabels, ['col1', 'col2', '3'])
        assert_equal([str, str, str], [type(x) for x in dt.clabels])

    def test_wrong_data_dimension(self):
        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels[:-1], self.clabels)
        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels, self.clabels[:-1])
        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels[:-1], self.clabels[:-1])

    def test_dtype(self):
        for dtype in (float, int, np.float32, np.bool, np.complex,
                      'float', 'int', 'float32', 'bool', 'complex'):
            dt = DataTable(self.data, self.rlabels, self.clabels, dtype=dtype)
            assert_equal(np.dtype(dtype), dt.dtype)

            dt = DataTable(self.data, self.rlabels, self.clabels,
                           dtype=np.dtype(dtype))
            assert_equal(np.dtype(dtype), dt.dtype)

        # Only builtin dtype accepted (no record-array!)
        assert_raises(PyCGHException, DataTable, self.data,
                      dtype=np.dtype([('f0', float)]))

    def test_dtypestr(self):
        dt = DataTable(self.data, dtype=str)
        assert_equal(np.dtype('S1'), dt.dtype)

    def test_reading(self):
        dt = DataTable(self.data)
        data = np.array(self.data)

        # rows
        for r in xrange(dt.rnum):
            assert_equal(data[r], dt[r])

        # columns
        for c in xrange(dt.cnum):
            assert_equal(data[:,c], dt[:,c])

        # both
        assert_equal(1, dt[0,0])
        assert_equal(6, dt[1,2])

        # slicing
        assert_equal(data, dt[0:2])
        assert_equal(data[0:1], dt[0:1])
        assert_equal(data[:,0:1], dt[:,0:1])
        assert_equal(data[:,0:2], dt[:,0:2])
        assert_equal(data[:,0:3], dt[:,0:3])
        assert_equal(data[0:1,0:2], dt[0:1,0:2])

    def test_label_reading(self):
        dt = DataTable(self.data)

        data = np.array(self.data)
        # Remember that:
        #      c0  c1  c2
        # r0    1   2   3
        # r1    4   5   6

        tests = ( # pairs 'index on ndarray', 'index on DataTable'
            # Row access
            (0, 0),
            (1, 1),
            (0, 'r0'),
            (1, 'r1'),
            (slice(0, 1, None), slice(0, 1, None)),
            (slice(0, None, None), slice(0, None, None)),
            (slice(0, 1, None), slice('r0', 'r1', None)),
            (slice(0, None, None), slice('r0', None, None)),

            # Column access
            ((0, slice(0, None, None)), ('r0', slice('c0', None, None))),
            ((0, slice(0, 2, None)), (0, slice('c0', 'c2', None))),
            ((0, slice(0, None, 2)), (0, slice('c0', None, 2))),

            # Combined access
            ( (slice(0, None, None), slice(0, None, None)),
              (slice('r0', None, None), slice('c0', None, None))),
            ( (slice(0, 1, None), slice(0, 2, None)),
              (slice('r0', 'r1', None), slice('c0', 'c2', None))),
            ( (slice(0, None, 2), slice(0, None, 2)),
              (slice('r0', None, 2), slice('c0', None, 2))),

            # List access
            ([0, 1], [0, 1]),
            ([0, 1], ['r0', 'r1']),
            (([0, 1, 1], [0, 1, 2]), (['r0', 'r1', 'r1'], ['c0', 'c1', 'c2'])),
        )

        for nd_test, dt_test in tests:
            assert_equal(data[nd_test], dt[dt_test])

    def test_index_of(self):
        dt = DataTable(self.data)

        assert_equal(0, dt.rindex('r0'))
        assert_equal(1, dt.rindex('r1'))
        assert_raises(PyCGHException, dt.rindex, 'r2')

        assert_equal(0, dt.cindex('c0'))
        assert_equal(1, dt.cindex('c1'))
        assert_raises(PyCGHException, dt.cindex, 'c3')

    def test_ndarray_conversion(self):
        dt = DataTable(self.data)
        a = np.asarray(dt)

        assert_equal(dt, a)
        a[0, 0] = 10
        assert_equal(dt, a)

    def test_load(self):
        dt_text = DataTable.load(StringIO(DATA_TEXT), delimiter=',')
        dt_data = DataTable(self.data)

        assert_equal(np.asarray(dt_text), np.asarray(dt_data))
        assert_equal(dt_text.rnum, dt_data.rnum)
        assert_equal(dt_text.cnum, dt_data.cnum)

    def test_save(self):
        dt = DataTable.load(StringIO(DATA_TEXT), delimiter=',')

        # Saving
        out = StringIO()
        dt.save(out)

        # Reading saved file
        dt_rel = DataTable.load(StringIO(out.getvalue()), delimiter=',')

        assert_equal(len(dt), len(dt_rel))
        assert_equal(dt.clabels, dt_rel.clabels)
        assert_equal(dt.rlabels, dt_rel.rlabels)
        assert_equal(np.asarray(dt), np.asarray(dt_rel))
