from collections import defaultdict

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2 import rinterface as ri

import rpy2.rlike.container as rlc

from ..datatypes.arraycgh import ArrayCGH

### CGHnormaliter objects -----------------------------------------------------
importr('CGHnormaliter')
CGHnormaliter = robjects.r['CGHnormaliter']
copynumber = robjects.r['copynumber']
calls = robjects.r['calls']
segmented = robjects.r['segmented']


def average_duplication(acgh, avg_funct=np.mean):
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


def cghnormaliter(acgh, nchrom=None, cellularity=1.0):
    # Averaging locations (id duplications) -> shrinked new data
    acgh = average_duplication(acgh)

    if nchrom is None:
        nchrom = len(np.unique(acgh['chromosome']))

    # Building an R DataFrame
    od = rlc.OrdDict([('CloneID', robjects.StrVector(acgh['id'])),
                      ('Chromosome', robjects.IntVector(acgh['chromosome'])),
                      ('Start', robjects.IntVector(acgh['start_base'])),
                      ('End', robjects.IntVector(acgh['end_base'])),
                      ('Test', robjects.FloatVector(acgh['test_signal'])),
                      ('Ref', robjects.FloatVector(acgh['reference_signal'])),
                      ('FakeTest', robjects.FloatVector(acgh['test_signal'])),
                      ('FakeRef', robjects.FloatVector(acgh['reference_signal'])),
                      ])
    df = robjects.DataFrame(od)
    result = CGHnormaliter(df, nchrom=nchrom, cellularity=cellularity,
                           plot_MA=False)

    # Attaching results
    acgh['cghnormaliter_ratio'] = np.asanyarray(copynumber(result))[:,0]
    acgh['cghnormaliter_call'] = np.asanyarray(calls(result), dtype=int)[:,0]
    acgh['cghnormaliter_segment'] = np.asanyarray(segmented(result))[:,0]
    return acgh
