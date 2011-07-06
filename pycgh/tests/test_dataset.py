#-*- coding: utf-8 -*-

# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

from nose.tools import *

import numpy as np
from pycgh.datatypes import Dataset
from pycgh import PyCGHException

class TestDataset(object):

    def test_creation(self):
        names = ('value1', 'value2')

        for t in [float, int, 'S10']:
            ds = Dataset(names, types=t)
            yield self.check_dataset, ds, names, (t,)*len(names)

    def test_creation_with_type(self):
        names = ('value1', 'value2')

        for t in [(float, float), (float, int), ('S10', int)]:
            ds = Dataset(names, types=t)
            yield self.check_dataset, ds, names, t

    def check_dataset(self, ds, names, types):
        assert_equal(0, len(ds))
        assert_equals(2, ds.dim)
        assert_equals(names, ds.names)
        assert_equals(types, ds.types)

    def test_creation_errors(self):
        assert_raises(PyCGHException, Dataset,
                                      ['value1', 'value2'],
                                      [float])
        assert_raises(PyCGHException, Dataset,
                                      ['value1', 'value2'],
                                      [float, float, float])
        assert_raises(PyCGHException, Dataset,
                                      ['value1', 'value2'],
                                      ['S10'])

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

    def test_auto_labels(self):
        ds = Dataset(['value1', 'value2'])
        ds.add(None, [1, 1])
        ds.add(None, [2, 2])

        assert_equals(('sample1', 'sample2'), ds.labels)

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
