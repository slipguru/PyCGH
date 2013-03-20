from collections import defaultdict

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2 import rinterface as ri

import rpy2.rlike.container as rlc

### CGHnormaliter objects -----------------------------------------------------
importr('CGHnormaliter')
CGHnormaliter = robjects.r['CGHnormaliter']
copynumber = robjects.r['copynumber']
calls = robjects.r['calls']
featureNames = robjects.r['featureNames']

def cghnormaliter(acgh, nchrom=None, cellularity=1.0, return_calls=False):

    # Filtered data extraction
    filtered = acgh.F[['id', 'chromosome', 'start_base', 'end_base',
                       'test_signal', 'reference_signal']]

    # Averaging locations
    if len(filtered['id']) > len(np.unique(filtered['id'])):
        id_indexes = defaultdict(list)

        for i, id in enumerate(filtered['id']):
            id_indexes[id].append(i)

        # First in order are kept
        selection = sorted([id_indexes[x][0] for x in id_indexes])
        for i in selection:
            probe_id = filtered['id'][i]
            idxs = id_indexes[probe_id]
            filtered[i]['test_signal'] = np.mean(filtered['test_signal'][idxs])
            filtered[i]['reference_signal'] = np.mean(filtered['reference_signal'][idxs])

        filtered = filtered[selection]
    else:
        selection = None

    if nchrom is None:
        nchrom = len(np.unique(filtered['chromosome']))

    # Building an R DataFrame
    od = rlc.OrdDict([('CloneID', robjects.StrVector(filtered['id'])),
                      ('Chromosome', robjects.IntVector(filtered['chromosome'])),
                      ('Start', robjects.IntVector(filtered['start_base'])),
                      ('End', robjects.IntVector(filtered['end_base'])),
                      ('Test', robjects.FloatVector(filtered['test_signal'])),
                      ('Ref', robjects.FloatVector(filtered['reference_signal'])),
                      ('FakeTest', robjects.FloatVector(filtered['test_signal'])),
                      ('FakeRef', robjects.FloatVector(filtered['reference_signal'])),
                      ])
    df = robjects.DataFrame(od)
    result = CGHnormaliter(df, nchrom=nchrom, cellularity=cellularity,
                           plot_MA=False)

    if return_calls:
        if selection:
            return selection, np.array(copynumber(result)), np.array(calls(result))
        return np.array(copynumber(result)), np.array(calls(result))
    else:
        if selection:
            return selection, np.array(copynumber(result))
        return np.array(copynumber(result))
