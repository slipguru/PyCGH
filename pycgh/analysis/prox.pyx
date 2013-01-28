import numpy as np
cimport numpy as np
np.import_array()

cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t

@cython.boundscheck(False)
def prox_squared_l1(np.ndarray[DTYPE_t, ndim=1] x not None,
                    np.float t not None):
    assert x.dtype == DTYPE
    
    cdef int n = x.shape[0]
    cdef int l
    cdef np.ndarray[DTYPE_t, ndim=1] abs_x = np.abs(x)
    cdef np.ndarray[DTYPE_t, ndim=1] x_s = np.sort(abs_x)[::-1]

    cdef np.ndarray[DTYPE_t, ndim=1] right_cond = x_s.cumsum()
    cdef np.ndarray[DTYPE_t, ndim=1] left_cond = 1./(2. * t) + np.arange(1, n+1)

    for l in range(n-1):
        if ( (left_cond[l] * x_s[l+1] <= right_cond[l]) and
             (left_cond[l] * x_s[l] > right_cond[l]) ):
            break

    cdef DTYPE_t rho = right_cond[l] / left_cond[l] # last computed conditions
    return np.sign(x) * np.clip(abs_x - rho, a_min=0.0, a_max=np.inf)