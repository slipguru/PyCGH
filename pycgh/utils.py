import itertools as it

import numpy as np

from rpy2 import robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri

from collections import defaultdict

def average_duplication(acgh, avg_funct=np.mean):
    """
    Returns a data array where every probe appears only once:
    for probes which are duplicated in the original array
    the value of the resulting probe is obtained averaging
    the original values.
    
    Parameters
    ----------
    
        acgh : :class:`pycgh.datatypes.ArrayCGH`
            The array CGH object to be modified.
            
        avg_funct : func, optional (default: numpy.mean)
            The function used to average values.
    """
    
    # Filtered data copy
    data = acgh.F[list(ArrayCGH.COL_NAMES)]

    # Averaging locations (id duplications)
    if len(data['id']) > len(np.unique(data['id'])):

        # Grouping duplications
        id_indexes = defaultdict(list)
        for i, id in enumerate(data['id']):
            id_indexes[id].append(i)

        # First position (in order) are kept
        selection = sorted([id_indexes[x][0] for x in id_indexes])
        for i in selection:
            probe_id = data['id'][i]
            idxs = id_indexes[probe_id]
            data[i]['test_signal'] = avg_funct(data['test_signal'][idxs])
            data[i]['reference_signal'] = avg_funct(data['reference_signal'][idxs])

        data = data[selection]

    return ArrayCGH(*(data[x] for x in ArrayCGH.COL_NAMES)) #sorted

def lowess(x, y, **kwargs):
    """
    Wrapper for the R implementation of the lowess algorithm
    """
    rlowess = ro.r['lowess']
    return rlowess(x, y, **kwargs)

def loess(values, x, y=None, **kwargs):
    """
    Wrapper for the R implementation of the loess algorithm
    """
    loess = ro.r['loess']
    predict = ro.r['predict']

    if y is None:
        fmla = ro.Formula('v ~ x')
        env = {'v': values, 'x': x}
        df = ro.DataFrame({'x': x})
    else:
        fmla = ro.Formula('v ~ x*y')
        env = {'v': values, 'x': x, 'y': y}
        df = ro.DataFrame({'x': x, 'y': y})

    for k in env:
        fmla.environment[k] = env[k]

    out = loess(fmla, **kwargs)
    return np.asarray(predict(out, df))

def array_trend(values, col, row):
    return loess(values, col, row, span=0.03, degree=1,
                 normalize=True, family='gaussian', iterations=3)

def probes_average(probes_id, probes_values, avg_function=np.mean):
    summary = dict()
    for id, v in it.izip(probes_id, probes_values):
        summary.setdefault(id, []).append(v)

    return dict((id, avg_function(summary[id])) for id in summary)

# IO utils --------------------------------------------------------------------
def _file_handle(file_ref, mode='r'):
    """
    Simply returns the file handler to `file_ref`, which can be either a string or a file handler. If the file to be opened is in a .gz format, it will be uncompressed.
    
    Parameters
    ----------
    
    file_ref : file or str
        Either a file handler or the path to the file which must be opened.
    
    mode : str, optional (default: 'r')
        The mode the file should be opened (read, write or append). If file_ref is already a file handler, thie parameter is ignored.
    
    Returns
    -------
    
    fh : file
        The file handler.
    """
    if not mode in 'rw':
        raise ValueError("mode must be 'r' or 'w'")

    def _is_string_like(obj):
        try:
            obj + ''
        except (TypeError, ValueError):
            return False
        return True

    try:
        if _is_string_like(file_ref):
            if file_ref.endswith('.gz'):
                import gzip
                fh = gzip.open(file_ref, mode='%sb' % mode)
            else:
                if mode == 'r':
                    fh = open(file_ref, 'U')
                else:
                    fh = open(file_ref, 'w')
        else:
            fh = file_ref
    except TypeError:
        raise ValueError('input file must be a path or file handle')

    return fh

# Splitting algoritm ----------------------------------------------------------
#from .datatypes.arraycgh import ArrayCGH

### WORKAROUND
ACGH_MISSING_INT = -1

def location_normalization(chr, start, end):
    # in some files the range is swapped :-/
    if start > end:
        start, end = end, start

    chr = chr.split('_', 1)[0].replace('chr', '') # from chrXX_bla_bla to XX
    if chr == 'X':
        return 23, start, end, False
    elif chr=='Y':
        return 24, start, end, False
    else:
        try:
            chr = int(chr)
            if chr <= 0 or chr > 24:
                # return (ArrayCGH.MISSING_INT,
                #         ArrayCGH.MISSING_INT,
                #         ArrayCGH.MISSING_INT, True)
            
                return (ACGH_MISSING_INT,
                        ACGH_MISSING_INT,
                        ACGH_MISSING_INT, True)
            
            return chr, start, end, False
        except ValueError: # unplaceable probe eg chrUn
            # return (ArrayCGH.MISSING_INT,
            #         ArrayCGH.MISSING_INT,
            #         ArrayCGH.MISSING_INT, True)
        
            return (ACGH_MISSING_INT,
                    ACGH_MISSING_INT,
                    ACGH_MISSING_INT, True)

def split_location(location):
    try:
        chr, interval = location.split(':')
        start, end = (int(x) for x in interval.split('-'))
    except ValueError: # unmapped/control probe
        # return (ArrayCGH.MISSING_INT,
        #         ArrayCGH.MISSING_INT,
        #         ArrayCGH.MISSING_INT, True)
        
        return (ACGH_MISSING_INT,
                ACGH_MISSING_INT,
                ACGH_MISSING_INT, True)
    
    return location_normalization(chr, start, end)