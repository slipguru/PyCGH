#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

from nose.tools import *

import numpy as np
from pycgh.datatypes import DataTable, Dataset
from pycgh import PyCGHException


class TestDataTable(object):
    def setup(self):
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.rlabels = ['r0', 'r1']
        self.clabels = ['c0', 'c1', 'c2']

    def test_creation(self):
        dt = DataTable(self.data, self.rlabels, self.clabels)
        assert_equal(2, dt.nrow)
        assert_equal(3, dt.ncol)
        assert_equal(2, len(dt))

        assert_equals(self.rlabels, dt.rlabels)
        assert_equals(self.clabels, dt.clabels)

    def test_labels(self):
        dt = DataTable(self.data)
        assert_equals(self.rlabels, dt.rlabels)
        assert_equals(self.clabels, dt.clabels)

        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels[:-1], self.clabels)
        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels, self.clabels[:-1])
        assert_raises(PyCGHException, DataTable, self.data,
                      self.rlabels[:-1], self.clabels[:-1])

    def test_dtype(self):
        for dtype in (float, int, np.float32, np.bool, np.complex):
            dt = DataTable(self.data, self.rlabels, self.clabels, dtype=dtype)
            assert_equals(dtype, dt.dtype)

            dt = DataTable(self.data, self.rlabels, self.clabels,
                           dtype=np.dtype(dtype))
            assert_equals(dtype, dt.dtype)

        # Only builtin dtype accepted
        assert_raises(PyCGHException, DataTable, self.data,
                      dtype=np.dtype([('f0', float)]))

    def test_dtypestr(self):
        dt = DataTable(self.data, dtype=str)
        assert_equals(np.asarray([1], dtype=str).dtype, dt.dtype)

    def test_reading(self):
        dt = DataTable(self.data)
        data = np.array(self.data)

        # rows
        assert_true(np.all(data[0] == dt[0]))
        assert_true(np.all(data[1] == dt[1]))

        #columns
        assert_true(np.all(data[:,0] == dt[:,0]))
        assert_true(np.all(data[:,1] == dt[:,1]))
        assert_true(np.all(data[:,2] == dt[:,2]))

        #both
        assert_true(np.all([1] == dt[0,0]))
        assert_true(np.all([6] == dt[1,2]))

        #slicing
        assert_true(np.all(data == dt[0:2]))
        assert_true(np.all(data[0:1] == dt[0:1]))
        assert_true(np.all(data[:,0:1] == dt[:,0:1]))
        assert_true(np.all(data[:,0:2] == dt[:,0:2]))
        assert_true(np.all(data[:,0:3] == dt[:,0:3]))
        assert_true(np.all(data[0:1,0:2] == dt[0:1,0:2]))

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
            assert_true(np.all(data[nd_test] == dt[dt_test]))




class TestDataset(object):


    def test_add_samples(self):
        ds = Dataset(['value1', 'value2'])

        # Adding
        ds.add('sample1', [1, 2])
        assert_equal(1, len(ds))

        # Reading
        assert_equal(1, ds['sample1']['value1'])
        assert_equal(2, ds['sample1']['value2'])
        assert_equal((1, 2), tuple(ds['sample1']))

        # Adding again
        ds.add('sample2', [2, 2])
        assert_equal(2, len(ds))

        # Reading again
        assert_equal(2, ds['sample2']['value1'])
        assert_equal(2, ds['sample2']['value2'])
        assert_equal((2, 2), tuple(ds['sample2']))

        # Failed adding
        assert_raises(ValueError, ds.add, 'sample1', [3, 3])

    def test_index_of(self):
        ds = Dataset(['value1', 'value2'])

        ds.add('sample1', [1, 1])
        ds.add('sample2', [2, 2])
        ds.add('sample10', [10, 10])
        ds.add(None, [3, 3]) # converted to sample4

        # Insert order
        assert_equals(0, ds.index_of('sample1'))
        assert_equals(1, ds.index_of('sample2'))
        assert_equals(2, ds.index_of('sample10'))
        assert_equals(3, ds.index_of('sample4'))

    def test_data_manually(self):

        for data in ([[1, 2], [3, 4]],
                     [(1, 2), (3, 4)],
                     ((1, 2), (3, 4)),
                     np.array([[1, 2], [3, 4]])):

            ds = Dataset(['value1', 'value2'],
                         data=data)
            yield self.check_given_data, ds

    def check_given_data(self, ds):
        assert_equals(2, len(ds))
        assert_equals(2, ds.dim)
        assert_equals(('value1', 'value2'), ds.names)
        assert_equals(('sample1', 'sample2'), ds.labels)

    def test_wrong_data(self):
        assert_raises(PyCGHException, Dataset,
                      ['value1', 'value2'], data=[[1, 2], [3, 4, 5]])
        assert_raises(PyCGHException, Dataset,
                      ['value1'], data=[[1, 2], [3, 4]])
        assert_raises(PyCGHException, Dataset,
                      ['value1'], data=[[1, 2], [3, 4, 5]])

    def test_array_conversion(self):
        ds = Dataset(['value1', 'value2'],
                     data=[[1, 2], [3, 4]])

        raw_data = ds.toarray()
        assert_equals((2,), raw_data.shape)

        print ds['value1']







class TestLabeledMatrix(object):
    #def _test_add_samples(self):


    def _test_delete_sample(self):
        lm = LabeledMatrix(['value1', 'value2'])
        for i in xrange(1, 11):
            lm.append('sample%d' % i, [i, i])

        for i in (2, 5, 9):
            assert_true(all([i, i] == lm['sample%d' % i]))

        idx = lm.index_of('sample10')
        assert_true(all([10, 10] == lm['sample10']))
        assert_true(all([10, 10] == lm[idx]))

        del lm['sample2']
        assert_raises(ValueError, lm.__getitem__, 'sample2')

        idx = lm.index_of('sample10')
        assert_true(all([10, 10] == lm['sample10']))
        assert_true(all([10, 10] == lm[idx]))

    def _test_change_sample(self):
        lm = LabeledMatrix(['value1', 'value2'])

        lm.append('sample1', [1, 2])
        lm['sample1'] = [5, 6]
        assert_equal(1, len(lm))

        del lm['sample1']
        assert_equal(0, len(lm))

        assert_raises(ValueError, lm.__delitem__, 'sample1')

    def _test_access_by_index(self):
        lm = LabeledMatrix(['value1', 'value2'])

        lm.append('sample1', [1, 1])
        lm.append('sample2', [2, 2])

        assert_true(all(lm[0] == lm['sample1']))
        assert_true(all(lm[1] == lm['sample2']))

        del lm['sample1']
        assert_true(all(lm[0] == lm['sample2']))



    def _test_load(self):
        lm = LabeledMatrix.load(os.path.join(PAR_DIR,
                                             'test_labeledmatrix.txt'))

        assert_equals(6, len(lm))
        assert_equals(['attr1', 'attr2'], sorted(lm.names))
        assert_true(all([10, 20] == lm['sample1']))
        assert_true(all([100, 200] == lm['sample2']))

    def _test_save(self):
        import tempfile
        out = os.path.join(tempfile.gettempdir(), 'lm.txt')

        lm = LabeledMatrix.load(os.path.join(PAR_DIR,
                                             'test_labeledmatrix.txt'))
        lm.save(out)
        lm2 = LabeledMatrix.load(out)

        assert_equals(len(lm), len(lm2))
        assert_equals(sorted(lm.names), sorted(lm2.names))
        assert_equals(sorted(lm.samples), sorted(lm2.samples))

        for sample in lm.samples:
            assert_true(all(lm[sample] == lm2[sample]))
