from cStringIO import StringIO

from nose.tools import *

from ..datatypes.cytobands import CytoStructure, ChromosomeStructure

CytoFileContent ="""\
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

#######
#def test_bands_structure():
#    for label in ('1', '1p', '1p1', '1p11', '1p11.2'):
#        cb = ChromosomeBand(label)
#        assert_equals(1, cb.chromosome)
#        assert_equals('1', cb.str_chromosome)
#        
#    assert_equals('X', ChromosomeBand('Xp11').str_chromosome)
#    assert_equals(23, ChromosomeBand('Xp11').chromosome)
#    
#    assert_raises(ValueError, ChromosomeBand, '1p112')
#    assert_raises(ValueError, ChromosomeBand, '1p11.')
#    assert_raises(ValueError, ChromosomeBand, 'p11')
#    assert_raises(ValueError, ChromosomeBand, '25p11')
#    assert_raises(ValueError, ChromosomeBand, '23p1.1')
#    
#def test_adding_bands():
#    chr1 = ChromosomeBand('1')
#    assert_equals(None, chr1.start_base)
#    assert_equals(None, chr1.end_base)
#    
#    chr1.update(ChromosomeBand('1p')) # arm
#    assert_equals(None, chr1.start_base)
#    assert_equals(None, chr1.end_base)
#    
#    chr1.update(ChromosomeBand('1p', 0, 100)) # arm
#    assert_equals(0, chr1.start_base)
#    assert_equals(100, chr1.end_base)
#
#### Forse e' inutilmente complicato!!

#####

def test_types():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    assert_equals(3, len(cs)) # number of chromosomes
    
    assert_equals(ChromosomeStructure, type(cs[1]))
    assert_equals(8, len(cs[1]))
    assert_equals(11, len(cs['Y']))
    assert_equals(list(cs[24]), list(cs['Y']))

def test_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    # Iteration is not dictionary like
    for chr, (int_chr, str_chr) in zip(cs, ((1, '1'), (23, 'X'), (24, 'Y'))):
        assert_equals(ChromosomeStructure, type(chr))
        assert_equals(int_chr, chr.chromosome)
        assert_equals(str_chr, chr.str_chromosome)

def test_chromosome_iteration():
    cs = CytoStructure(StringIO(CytoFileContent))

    # standard iteration
    assert_equals(['p12', 'p11.2', 'p11.1', 'q11',
                   'q12', 'q21.1', 'q21.2', 'q21.3'], list(cs[1]))
    
    # resolution based iteration
    assert_equals(list(cs[1]), list(cs[1].bands_iter()))
    assert_equals(['p1', 'q1', 'q2'], list(cs[1].bands_iter(level=1)))
    assert_equals(['p12', 'p11', 'q11',
                   'q12', 'q21'], list(cs[1].bands_iter(level=2)))
    
    assert_equals(list(iter(cs[1])), list(cs[1].bands_iter(level=3))) # full
    assert_equals(list(iter(cs[1])), list(cs[1].bands_iter(level=4))) # out
       
def test_slicing():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    assert_equals(['p12', 'p11.2'], cs[1][120700000:121100000])
    assert_equals(['p11.2'], cs[1][120700001:121100000])
    assert_equals(['q11.21'], cs['Y'][13500000:14000000])
    
def test_slicing_notucsc():
    NotUCSCContent_ = ("chr1 117600001 120700000 p12   gpos50\n"
                       "chr1 120700001 121100000 p11.2 gneg\n"
                       "chr1 121100001 124300000 p11.1 acen")
    
    cs = CytoStructure(StringIO(NotUCSCContent_), format=None)
    assert_equals(['p12', 'p11.2'], cs[1][120700000:121100000])
    assert_equals(['p11.2'], cs[1][120700001:121100000])

def test_slicing_defaults():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    assert_equals(7, len(cs[1][120700001:]))
    assert_equals(7, len(cs[1][:149600000]))
    
    assert_equals(8, len(cs[1][:])) #ALL
    assert_equals(8, len(cs['1'])) #ALL
    assert_equals(11, len(cs['Y'][:])) #ALL

def test_positions():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    assert_equals((117600001, 153300000), cs[1].position())
    
    # Descresing resolution
    assert_equals((14300001, 19000000), cs['Y'].position('q11.221'))
    assert_equals((14300001, 25400000), cs['Y'].position('q11.22'))
    assert_equals((12500001, 27200000), cs['Y'].position('q11.2'))
    assert_equals((11300001, 27200000), cs['Y'].position('q11'))
    assert_equals((11300001, 57772954), cs['Y'].position('q1'))
    assert_equals((11300001, 57772954), cs['Y'].position('q'))
    
    # Default
    assert_equals(cs[1].position(), cs[1].position(''))
    
    # Errors
    assert_raises(ValueError, cs['Y'].position, 'f')
    assert_raises(ValueError, cs['Y'].position, '2')
    assert_raises(ValueError, cs['Y'].position, 2)
    assert_raises(ValueError, cs['Y'].position, 'p7')
    assert_raises(ValueError, cs['Y'].position, 'q1.')
    assert_raises(ValueError, cs['Y'].position, 'q11.')
    assert_raises(ValueError, cs['Y'].position, 'q2')

def test_expand_band():
    cs = CytoStructure(StringIO(CytoFileContent))
    
    # From pter to centroid
    assert_equals(['p11.32', 'p11.31', 'p11.2', 'p11.1'], cs['Y'].expand('p1'))
    
    # From centroid to qter
    assert_equals(['q11.1', 'q11.21', 'q11.221', 'q11.222', 
                   'q11.223', 'q11.23', 'q12'], cs['Y'].expand('q1'))
    
    # With position
    bands_table = cs['Y'].expand('p1', with_position=True)
    assert_equals(4, len(bands_table))
    assert_equals(tuple(cs['Y'].expand('p1')), zip(*bands_table)[0])

    # Errors
    assert_raises(ValueError, cs[1].expand, 'p36.')
    assert_raises(ValueError, cs[1].expand, 'z36')
    
def test_gstrand():
    chr1 = CytoStructure(StringIO(CytoFileContent))[1]
    
    map = {'p12': 'gpos50', 'p11.2': 'gneg', 'p11.1': 'acen', 'q21.3': 'gneg',
           'q11': 'acen', 'q12': 'gvar', 'q21.1': 'gneg', 'q21.2': 'gpos50'}
    
    for band in chr1:
        assert_equals(map[band], chr1.gstrand(band))
        
       
        