from libc.math cimport fabs
cimport numpy as np
import numpy as np

cimport cython
np.import_array()

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

cdef inline double fmax(double x, double y):
    if x > y:
        return x
    return y

cdef inline double fsign(double f):
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def prox_squared_l1(np.ndarray[DOUBLE, ndim=1] x, DOUBLE t):
    cdef unsigned int n = x.shape[0]
    cdef unsigned int l = 0

    cdef double cum = 0.0
    cdef double k = (1./(2. * t)) - n+2
    cdef double rho = 0.0

    cdef np.ndarray[DOUBLE, ndim=1] x_s = np.abs(x)
    x_s.sort()

    for l in range(n-1, 0, -1): # reversed order
        cum += x_s[l]
        if ((k+l) * x_s[l-1] <= cum) and ((k+l) * x_s[l] > cum):
            break

    rho = cum / ((k+l) * x_s[l]) # last computed conditions
    for l in range(n):
        x_s[l] = fsign(x[l]) * fmax(fabs(x[l]) - rho, 0.0)
    return x_s
