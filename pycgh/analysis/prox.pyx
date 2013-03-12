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
    cdef double k = (1./(2. * t)) + n
    cdef double rho = 0.0

    cdef np.ndarray[DOUBLE, ndim=1] x_s = np.abs(x)
    if n == 1:
        l = 0
        cum = x_s[l]
    else:
        x_s.sort()
    
        for l in range(n-1, 0, -1): # reversed order
            cum += x_s[l]
            if ((k-l) * x_s[l-1] <= cum) and ((k-l) * x_s[l] > cum):
                break

    rho = 2.*t*cum / (1 + 2*t*(n-l))
    for l in range(n):
        x_s[l] = fsign(x[l]) * fmax(fabs(x[l]) - rho, 0.0)
    return x_s

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def prox_squared_l1_bycol(np.ndarray[DOUBLE, ndim=2] X,
                          np.ndarray[DOUBLE, ndim=2] Y,
                          DOUBLE t):
    cdef unsigned int num_row = X.shape[0]
    cdef unsigned int num_col = X.shape[1]
    cdef unsigned int l = 0
    cdef unsigned int c = 0

    cdef double cum = 0.0
    cdef double rho = 0.0
    cdef double k = (1./(2. * t)) + num_row
       
    Y[:,:] = np.abs(X)
    for c in range(num_col):
        if num_row == 1:
            l = 0
            cum = Y[0, c]
        else:
            Y[:,c].sort()
    
            cum = 0.0
            for l in range(num_row-1, 0, -1): # reversed order (exit with l=0)
                cum += Y[l,c]
                if ((k-l) * Y[l-1,c] <= cum) and ((k-l) * Y[l,c] > cum):
                    break

        rho = 2.*t*cum / (1 + 2*t*(num_row-l))
        for l in range(num_row):
            Y[l,c] = fsign(X[l,c]) * fmax(fabs(X[l,c]) - rho, 0.0)
