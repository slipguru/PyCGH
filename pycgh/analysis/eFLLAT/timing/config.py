# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Configuration file imported at run-time """

import sys
import os

                                 # select the python module:
ALGORITHMS = ['fllat_nowak',     # 0
              'fllat',           # 1
              'l12_B',           # 2
              'tvw_B',           # 3
              'boxpos_Theta',    # 4
              'cghdl']           # 5
ALGORITHM = sys.argv[2]          # <--
if not ALGORITHM in ALGORITHMS:
    raise Exception('Algorithm %s not recognized!' % ALGORITHM)

# Parameters (Valid for all versions)
params = {
    'J_range': [3, 5, 7],
    'lambda_range': [1e-4, 1e-3, 1e-2, 1e-1],
    'mu_range': [1e0, 5],# 1e1],
    'maxK': 1e2, #external
    'maxN': 5e5, #internal
    'eps': 0.0,
    'eps_gap': 1e-3, # if 'auto': C/k
    'step': 1.0, # or 'decr', 'incr'
    'eps_jumps': 0.0
}

OUTPUT = os.path.expanduser('~')

# Parameters specific for the give module (model/algorithm)
if ALGORITHM == 'fllat_nowak':
    params['tvw'] = None                # **Ignored**
    params['theta_bound'] = None        # **Ignored**
    params['tau_range'] = [0.0]         # **Ignored**
if ALGORITHM == 'fllat':
    params['tvw'] = None                # **Ignored**
    params['theta_bound'] = 1.0         # l2 norm bound (nowak uses 1.0)
    params['tau_range'] = [0.0]         # **Ignored**
elif ALGORITHM == 'l12_B':
    params['tvw'] = None                # **Ignored**
    params['theta_bound'] = 1.0         # l2 norm bound (nowak uses 1.0)
    params['tau_range'] = [0.0]         # **Ignored**
elif ALGORITHM == 'tvw_B':
    params['tvw'] = 'auto'              # based on dataset infos
    params['theta_bound'] = 1.0         # l2 norm bound (nowak uses 1.0)
    params['tau_range'] = [0.0]         # **Ignored**
elif ALGORITHM == 'boxpos_Theta':
    params['tvw'] = 'auto'              # based on dataset infos
    params['theta_bound'] = 1.0         # BOX bound
    params['tau_range'] = [0.0]         # **Ignored**
elif ALGORITHM == 'cghdl':
    params['tvw'] = 'auto'              # based on dataset infos
    params['theta_bound'] = 1.0         # BOX bound
    params['tau_range'] = [1e-3, 1e-2]
                           #1e-1, 1e0]   # l12 theta range
