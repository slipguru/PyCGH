#=======================================
#from mlabwrap import mlab
#mlab.addpath('chen/')
#mlab.addpath('chen/stprtool/')
#mlab.stprpath('chen/stprtool/')
#data = mlab.ReadAgilentResult(FILE_PATH)
#aCGHdata, model = mlab.aCGHNormalization3(data, 1, (0 if swap else 1), 1, nout=2)
#plt.scatter(np.arange(len(aCGHdata.nRatioAdj)),
#            aCGHdata.nRatioAdj.flatten(),
#            c=aCGHdata.nRatioAdj.flatten(),
#            cmap=plt.get_cmap('jet'),
#            vmin=-1, vmax=1, s=8, edgecolors='none')
#plt.show()
#exit()
#=======================================


## Manual lowess normalization
#if show_lowess:
#    rlowess = robjects.r['lowess'] # Lowess from R
#    lowess_curve = rlowess(A, M, f=2./3., iter=3) #std params 2/3
#    x = np.asarray(lowess_curve[0])
#    y = np.asarray(lowess_curve[1])
#    plt.plot(x, y, 'r-', lw=2, label='%s Lowess Curve' % label) # Lowess Line
#
#    # M must be updated sorted
#    assert np.allclose(x, np.sort(A))
#    sorted_idxs = np.argsort(A) # as 'x' in lowess_curve
#    M_norm = np.empty_like(M)
#    M_norm[sorted_idxs] = M[sorted_idxs] - y
#
#    plt.scatter(A, M_norm, label='Lowess-Normalized %s' % label,
#                s=4, edgecolors='none', c='c')
#    plt.axhline(np.median(M_norm), lw=2, c='c',
#                label='Lowess-Normalized %s Median' % label)
#
#    lowess_curve = rlowess(A, M_norm, f=2./3, iter=3) #std params
#    x = np.asarray(lowess_curve[0])
#    y = np.asarray(lowess_curve[1])
#    plt.plot(x, y, 'b-', lw=2,
#             label='Lowess-Normalized %s Lowess Curve' % label)
#
#    return M_norm
