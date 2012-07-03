from cStringIO import StringIO

from nose.tools import *

from ..datatypes.cytobands import CytoStructure
from ..datatypes.cytobands import ChromosomeStructure, ChromosomeBand

CytoFileContent = """\
chr1	117600000	120700000	p12	gpos50
chr1	120700000	121100000	p11.2	gneg
chr1	121100000	124300000	p11.1	acen
chr1	149600000	153300000	q21.3	gneg
chr1	124300000	128000000	q11	acen
chr1	128000000	142400000	q12	gvar
chr1	142400000	148000000	q21.1	gneg
chr1	148000000	149600000	q21.2	gpos50
chrY	0	1700000	p11.32	gneg
chrY	1700000	3300000	p11.31	gpos50
chrY	3300000	11200000	p11.2	gneg
chrY	11200000	11300000	p11.1	acen
chrY	11300000	12500000	q11.1	acen
chrY	12500000	14300000	q11.21	gneg
chrY	14300000	19000000	q11.221	gpos50
chrY	19000000	21300000	q11.222	gneg
chrY	21300000	25400000	q11.223	gpos50
chrY	25400000	27200000	q11.23	gneg
chrY	27200000	57772954	q12	gvar
chrX	37500000	42300000	p11.4	gneg
chrX	42300000	47300000	p11.3	gpos75
chrX	47300000	49700000	p11.23	gneg
chrX	49700000	54700000	p11.22	gpos25
chrX	54700000	56600000	p11.21	gneg
chrX	56600000	59500000	p11.1	acen
chrX	59500000	65000000	q11.1	acen
chrX	65000000	65100000	q11.2	gneg
"""

def test_bands_ordering():
    cs = CytoStructure(StringIO(CytoFileContent))

    assert_true(cs[1].band('p12') < cs[1].band('p11.2'))
    assert_true(cs[1].band('p12') < cs[1].band('p11'))
    assert_true(cs[1].band('p') < cs[1].band('q'))
    assert_true(cs[1].band('q11') < cs[1].band('q2'))
    assert_true(cs[24].band('p') < cs['Y'].band('q'))

    # Equalities
    assert_true(cs[1].band('p12') == cs[1].band('p12'))
    assert_true(cs[1].band('p') == cs[1].band('p'))
    assert_true(cs[1].band('q') == cs[1].band('q'))
    assert_false(cs[1].band('p11.1') == cs['Y'].band('p11.1'))
    assert_true(cs[1].band('p11.1') != cs['Y'].band('p11.1'))

    # Overlapping bands
    assert_raises(RuntimeError, cs[1].band('p12').__cmp__, cs[1].band('p1'))
    assert_raises(RuntimeError, cs[1].band('p12').__cmp__, cs[1].band('p'))
    assert_raises(RuntimeError, cs[1].band('p').__cmp__, cs[1].band('p1'))

    # Different Chromosomes (natural ordering by chromosomes)
    assert_true(cs[1].band('p') < cs['Y'].band('q'))

def test_types():
    cs = CytoStructure(StringIO(CytoFileContent))

    assert_equals(3, len(cs)) # number of chromosomes

    assert_equals(ChromosomeStructure, type(cs[1]))
    assert_equals(8, len(cs[1]))
    assert_equals(11, len(cs['Y']))
    assert_equals(cs[24], cs['Y']) # ChromosomeStructure
    assert_equals(list(cs[24]), list(cs['Y'])) # list of ChromosomeBands

def test_get_band():
    cs = CytoStructure(StringIO(CytoFileContent))

    cb = ChromosomeBand(1, 'p12', 117600001, 120700000, 'gpos50')
    assert_equals(cb, cs[1].band('p12'))

    cb = ChromosomeBand(1, 'p1', 117600001, 124300000)
    assert_equals(cb, cs[1].band('p1'))

    # Default
    assert_equals(cs[1].band(), cs[1].band(''))

    # Errors
    assert_raises(ValueError, cs['Y'].band, 'f')
    assert_raises(ValueError, cs['Y'].band, '2')
    assert_raises(ValueError, cs['Y'].band, 2)
    assert_raises(ValueError, cs['Y'].band, 'p7')
    assert_raises(ValueError, cs['Y'].band, 'q1.')
    assert_raises(ValueError, cs['Y'].band, 'q11.')
    assert_raises(ValueError, cs['Y'].band, 'q2')

def test_get_band_values():
    cs = CytoStructure(StringIO(CytoFileContent))
    def pos(band): return band.start_base, band.end_base

    # Descresing resolution
    assert_equals((14300001, 19000000), pos(cs['Y'].band('q11.221')))
    assert_equals((14300001, 25400000), pos(cs['Y'].band('q11.22')))
    assert_equals((12500001, 27200000), pos(cs['Y'].band('q11.2')))
    assert_equals((11300001, 27200000), pos(cs['Y'].band('q11')))
    assert_equals((11300001, 57772954), pos(cs['Y'].band('q1')))
    assert_equals((11300001, 57772954), pos(cs['Y'].band('q')))

def test_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))

    # Iteration is not dictionary like
    for chr, (int_chr, str_chr) in zip(cs, ((1, '1'), (23, 'X'), (24, 'Y'))):
        assert_equals(ChromosomeStructure, type(chr))
        assert_equals(int_chr, chr.chromosome)
        assert_equals(str_chr, chr.str_chromosome)

def test_chromosome_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))

    # Reading from file data
    chr1_bands = [line for line in CytoFileContent.split('\n')
                                    if line.startswith('chr1')]
    data = (band.split()[1:] for band in chr1_bands)
    data = dict((d[2], {'sb': int(d[0]), 'eb': int(d[1]), 'gs': d[3]}) for d in data)

    for band in cs[1]:
        assert_equals(ChromosomeBand, type(band))
        assert_true(band.label in data)
        assert_equals(data[band.label]['sb'] + 1 , band.start_base)
        assert_equals(data[band.label]['eb'] , band.end_base)
        assert_equals(data[band.label]['gs'] , band.gstrand)

def test_chromosome_resolution_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))

    # resolution based iteration
    assert_equals(list(cs[1]), list(cs[1].bands_iter()))
    assert_equals([('p1', 117600001, 124300000),
                   ('q1', 124300001, 142400000),
                   ('q2', 142400001, 153300000)],
                  [(b.label, b.start_base, b.end_base)
                    for b in cs[1].bands_iter(level=1)])
    assert_equals([('p12', 117600001, 120700000),
                   ('p11', 120700001, 124300000),
                   ('q11', 124300001, 128000000),
                   ('q12', 128000001, 142400000),
                   ('q21', 142400001, 153300000)],
                  [(b.label, b.start_base, b.end_base)
                    for b in cs[1].bands_iter(level=2)])

    assert_equals(list(cs[1]), list(cs[1].bands_iter(level=3))) # full
    assert_equals(list(cs[1]), list(cs[1].bands_iter(level=4))) # out

    # gstrand test: initialized only for full resolution
    assert_equals(['gpos50', 'gneg', 'acen', 'acen',
                   'gvar', 'gneg', 'gpos50', 'gneg'],
                  [b.gstrand for b in cs[1].bands_iter()])
    assert_equals([None]*3, [b.gstrand for b in cs[1].bands_iter(level=1)])

def test_slicing():
    cs = CytoStructure(StringIO(CytoFileContent))

    assert_equals(['p12', 'p11.2'], [b.label for b in cs[1][120700000:121100000]])
    assert_equals(['p11.2'], [b.label for b in cs[1][120700001:121100000]])
    assert_equals(['q11.21'], [b.label for b in cs['Y'][13500000:14000000]])

def test_slicing_notucsc():
    NotUCSCContent_ = ("chr1 117600001 120700000 p12   gpos50\n"
                       "chr1 120700001 121100000 p11.2 gneg\n"
                       "chr1 121100001 124300000 p11.1 acen")

    cs = CytoStructure(StringIO(NotUCSCContent_), format=None)
    assert_equals(['p12', 'p11.2'], [b.label for b in cs[1][120700000:121100000]])
    assert_equals(['p11.2'], [b.label for b in cs[1][120700001:121100000]])

def test_slicing_defaults():
    cs = CytoStructure(StringIO(CytoFileContent))

    assert_equals(7, len(cs[1][120700001:]))
    assert_equals(7, len(cs[1][:149600000]))
    assert_raises(ValueError, cs[1].__getitem__, slice(0, None))

    assert_equals(8, len(cs[1][:])) #ALL
    assert_equals(8, len(cs['1'])) #ALL
    assert_equals(11, len(cs['Y'][:])) #ALL

def test_expand_band():
    cs = CytoStructure(StringIO(CytoFileContent))

    # From pter to centroid
    assert_equals(['p11.32', 'p11.31', 'p11.2', 'p11.1'],
                   [b.label for b in cs['Y'].band('p1').expand()])

    # From centroid to qter
    assert_equals(['q11.1', 'q11.21', 'q11.221', 'q11.222',
                   'q11.223', 'q11.23', 'q12'],
                  [b.label for b in cs['Y'].band('q1').expand()])

    assert_equals(['p12', 'p11.2', 'p11.1', 'q11',
                   'q12', 'q21.1', 'q21.2', 'q21.3'],
                  [b.label for b in cs[1].band().expand()])

    # Unexpandable band
    assert_equals(None, cs[1].band('p12').expand())

def test_expand_with_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))

    for band in cs[1].bands_iter(level=1):
        assert_equals(2, len(band.label))

        sub_bands = band.expand()
        assert_equals(sorted(sub_bands), list(sub_bands))
        for sub_band in sub_bands:
            assert_true(sub_band.label.startswith(band.label))

def test_gstrand():
    chr1 = CytoStructure(StringIO(CytoFileContent))[1]

    map = {'p12': 'gpos50', 'p11.2': 'gneg', 'p11.1': 'acen', 'q21.3': 'gneg',
           'q11': 'acen', 'q12': 'gvar', 'q21.1': 'gneg', 'q21.2': 'gpos50'}

    for band in chr1:
        assert_equals(map[band.label], band.gstrand)
        assert_equals(chr1.band(band.label).gstrand, band.gstrand)
