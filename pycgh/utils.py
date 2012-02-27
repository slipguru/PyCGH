import itertools as it

import numpy as np

from rpy2 import robjects
import rpy2.robjects.numpy2ri

def lowess(x, y, **kwargs):
    """ Lowess from R """
    rlowess = robjects.r['lowess']
    return rlowess(x, y, **kwargs)

def loess(values, x, y=None, **kwargs):
    """ Loess from R """
    loess = robjects.r['loess']
    predict = robjects.r['predict']

    if y is None:
        fmla = robjects.Formula('v ~ x')
        env = {'v': values, 'x': x}
        df = robjects.DataFrame({'x': x})
    else:
        fmla = robjects.Formula('v ~ x*y')
        env = {'v': values, 'x': x, 'y': y}
        df = robjects.DataFrame({'x': x, 'y': y})

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


#-----------
class ArrayCGHSynth(object):
    pass
