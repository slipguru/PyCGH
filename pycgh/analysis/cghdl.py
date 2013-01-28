# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD

import numpy as np

# CGH DL utils ----------------------------------------------------------------
## Projections -----
def pos_proj(x, r=np.inf):
    return np.clip(x, 0.0, r)

def ball_proj(x, r=1.0):
    y = np.array(x, dtype=float) # copy

    y_norm = np.linalg.norm(y)
    if y_norm > r:
        y *= (r / y_norm)

    return y

def symplex_proj(x):
    x_s = np.sort(x)[::-1]
    xtmp = (x_s.cumsum() - 1) / np.arange(1, len(x) + 1)
    k = np.where(xtmp < x_s)[0][-1]

    return np.maximum(x - np.max(xtmp[k], 0), 0)

class CGHDL(object):
    def __init__(self):
        pass
