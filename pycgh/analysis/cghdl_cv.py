import sys
import os
import itertools as it
import datetime as dt
import shutil
import logging

import numpy as np

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
robjects.conversion.py2ri = numpy2ri.numpy2ri

import pylab as pl
import matplotlib.colors

### Useful functions ----------------------------------------------------------
#from pycgh.analysis import cghDL
from cghDL import cghDL
from pycgh.datatypes.cytobands import CytoStructure
from pycgh.datatypes.dataset import DataTable
from pycgh.utils import probes_average

# we dont want to penalize slight variations within the same level
def atoms_jumps(B):
    L, J = B.shape
    jumps = 0
    for j in xrange(J):

        # How much different levels
        jumps += len(np.nonzero(np.unique(B[:,j]))[0])
    return jumps

#counts jumps bigger than eps tolerance
def approx_jumps(B, eps=1e-3):
    L, J = B.shape
    ajumps = 0
    BS=np.sort(B,axis=0)

    for j in xrange(J):
      db=np.diff(BS[:,j])
      ajumps += sum(db>eps)
    return ajumps

def plot_results(Y, B, Theta, J, alpha=None, lambda_=None, tau=None, mu=None,
                 jumpsB=None, fit=None, bp=None):
    fig = pl.figure()
    # Color boundaries
    vmin = -2.5 * Y.std()
    vmax = 2.5 * Y.std()
    norm = matplotlib.colors.normalize(vmax=vmax, vmin=vmin)

    if not bp is None:
        breaks = np.where(bp==0)[0]
    else:
        breaks = []

    if tau is None:
        pl.suptitle('$J=$%d - $\\alpha=$%6.10f - $\\lambda=$%6.10f ( $\\lambda_1=$%6.10f - $\\lambda_2=$%6.10f)' % (J, alpha, lambda_, alpha*lambda_, (1-alpha)*lambda_))
    else:
        pl.suptitle('$J=$%d - $\\lambda=$%6.10f - $\\mu=$%6.10f - $\\tau=$%6.10f' % (J, lambda_, mu, tau))

    pl.subplot(221)
    pl.imshow(Y, aspect='auto', interpolation='none', norm=norm)
    for b in breaks:
        pl.axhline(b+1, lw=.5, ls='-', c='k')
    pl.title('Y')

    pl.subplot(222)
    pl.imshow(np.dot(B, Theta), aspect='auto', interpolation='none', norm=norm)
    for b in breaks:
        pl.axhline(b+1, lw=.5, ls='-', c='k')
    if fit is None:
        pl.title('Y Estimated')
    else:
        pl.title('Y Est. (fit=%.4f)' % fit)

    pl.subplot(223)
    im_ref = pl.imshow(Theta, aspect='auto', interpolation='none', norm=norm)
    pl.title('Theta')

    pl.subplot(224)
    pl.imshow(B, aspect='auto', interpolation='none', norm=norm)
    for b in breaks:
        pl.axhline(b+1, lw=.5, ls='-', c='k')
    if jumpsB is None:
        pl.title('B')
    else:
        pl.title('B(cplx=%s)' % str(int(jumpsB)))

    # Make an axis for the colorbar on the right side
    pl.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(im_ref, cax=cax, extend='both')
### ---------------------------------------------------------------------------

if __name__ == '__main__':
    INPUT_DATA = sys.argv[1]
    OUTPUT_DIR = sys.argv[2]

    # Reading data, CGH Array Platform and Cytoband Structures
    ds = DataTable.load(INPUT_DATA, delimiter=',')
    dm = DataTable(np.loadtxt('agilentCgh4x44k.txt.gz', dtype=str)[:,1:],
                   clabels=['chr', 'start', 'end', 'id'])
    cs = CytoStructure('cytoBandIdeo.txt.gz')

    probes_starts = dict()
    for (p, c, s) in dm[:,['id', 'chr', 'start']]:
        try:
            c = c[3:]
            if c == 'X':
                c = 23
            elif c == 'Y':
                c = 24
            else:
                c = int(c)
            s = int(s)
            probes_starts[p] = (c, s)
        except ValueError:
            pass

    pairs = sorted([(probes_starts[probe], i) for i, probe in enumerate(ds.rlabels)])
    sort_idxs = [x[1] for x in pairs]

    data = ds[sort_idxs]
    chr_idxs = [i for i, p in enumerate(np.array(ds.rlabels)[sort_idxs]) if probes_starts[p][0] in (13, 15, 18, 21)]

    probes = np.array(ds.rlabels)[sort_idxs]
    probes = probes[chr_idxs]
    data = data[chr_idxs]

    bands = [str(cs[probes_starts[p][0]][probes_starts[p][1]:probes_starts[p][1]][0]).split(' ')[0] for p in probes]

    np.savez_compressed(os.path.join(OUTPUT_DIR, 'PreprocessedData.npz'),
                        data=data,
                        probes=probes,          # Row
                        bands=bands,
                        samples=ds.clabels)     # Columns

    # Breackpoints calculation (TV)
    chrl = np.asarray([band.split('p')[0] if 'p' in band else band.split('q')[0]
                       for band in bands])
    chrarmsl = np.asarray(['p' if 'p' in band else 'q' for band in bands])

    w = np.logical_and(chrl[1:] == chrl[:-1],
                       chrarmsl[1:] == chrarmsl[:-1])
    w = np.asarray(w, dtype=float)

    # Useful for plotting
    bp = np.asarray(chrl[1:] == chrl[:-1], dtype=float)

    Y = data
    L, S = Y.shape

    # PCA
    TMP = 1./np.sqrt(S-1) * Y.T
    TMP -= np.mean(TMP, axis=0)
    U, d, Vt = np.linalg.svd(TMP, full_matrices=False)
    V = np.array(Vt.T)
    #variance = d**2
    #print (np.cumsum(variance)/np.sum(variance))*100.
    #J_values = [22, 42, 55] # 50%, 70%, 80% varianza
    #J_values = [15] # more than 50% of variance
    J_values = [5] # sappiamo che sono 3 gruppi... andiamo larghi
    #exit()

    mu_values = np.logspace(0, 1, 2)      # B TV
    lambda_values = np.logspace(-3, -2, 2)  # B l1
    tau_values =  np.logspace(-2, -1, 2)     # T l1
    ###########################################################################

    ### Logging ---------------------------------------------------------------
    logger = logging.getLogger('cghDL')
    logger.setLevel(logging.DEBUG)
    logfilename = 'CGHDL_aCGH_%s.log' % dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    lfile = logging.FileHandler(os.path.join(OUTPUT_DIR, logfilename))
    lfile.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    lfile.setFormatter(formatter)
    logger.addHandler(lfile)
    ### -----------------------------------------------------------------------

    BICS = list()
    BIC_min = np.inf
    parameters = list(it.product(J_values, mu_values, lambda_values, tau_values))
    for J, mu, lambda_, tau in parameters:
        infomsg = ('    * Experiment with J = %2d, '
                   'mu = %6.10f, lambda = %6.10f,  '
                   'tau = %6.10f' % (J, mu, lambda_, tau))
        print infomsg
        logger.info(infomsg)

        result_cghdl = cghDL(Y, J=J, lambda_=lambda_, mu=mu,
                             tau=tau, tvw=w)#, init=V[:, :J])#, maxK=10, maxN=50)  # B init
        B = result_cghdl['B']
        Theta = result_cghdl['Theta']
        conv_iter = result_cghdl['conv']
        gap_phi = result_cghdl['gap_phi']
        gap_psi = result_cghdl['gap_psi']

        #jumpsB = atoms_jumps(B)
        jumpsB = approx_jumps(B)
        fit = np.sum((Y - np.dot(B, Theta))**2.) / (S*L)

        BIC = ( ( (S*L) * np.log(fit) )  +
                ( jumpsB * np.log(S*L)     ) ) # complexity
        BICS.append(BIC)

        infomsg = ('\t| conv at = %d' % conv_iter +
                   ', last gaps =(%.4f,%.4f)' %(gap_phi,gap_psi) +
                   ', BIC = %10.4f' % BIC)
        print infomsg
        logger.info(infomsg)

        plot_results(Y, B, Theta, J, mu=mu, lambda_=lambda_, tau=tau,
                     jumpsB=jumpsB, fit=fit, bp=bp)

        img_file = 'J%2d_M%6.10f_L%6.10f_T%6.10f_out.png' % (J, mu, lambda_, tau)
        infomsg = 'complete filename %s ' % os.path.join(OUTPUT_DIR, img_file)
        print infomsg
        logger.info(infomsg)
        pl.savefig(os.path.join(OUTPUT_DIR, img_file))

        pl.clf()
        pl.close()

        if BIC < BIC_min:
            np.savez_compressed(os.path.join(OUTPUT_DIR, 'Result.npz'),
                                Y=Y, B=B, Theta=Theta)
            BIC_min = BIC

    best = np.argmin(BICS)
    J, mu, lambda_, tau = parameters[best]
    infomsg = '    *** Best -> J = %2d, mu = %6.10f, lambda = %6.10f, tau = %6.10f, BIC min = %.3f (%.3f)' % (J, mu, lambda_, tau, BICS[best], BIC_min)
    print infomsg
    logger.info(infomsg)

    img_file = 'J%2d_M%6.10f_L%6.10f_T%6.10f_out.png' % (J, mu, lambda_, tau)
    shutil.copy2(os.path.join(OUTPUT_DIR, img_file), os.path.join(OUTPUT_DIR, 'Result.png'))

    BICresult = np.load(os.path.join(OUTPUT_DIR, 'Result.npz'))
    Y = BICresult['Y']
    B = BICresult['B']
    Theta = BICresult['Theta']
    infomsg = '    *** Saving best results (compressed .npz)...'
    print infomsg
    logger.info(infomsg)
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'Result.npz'),
                        Y=Y, B=B, Theta=Theta,
                        parameters=parameters,
                        BICS=BICS)
    # un po' ridondante...

    pl.figure()
    pl.plot(BICS, '.-')
    pl.savefig(os.path.join(OUTPUT_DIR, 'BICS.png'))

    pl.clf()
    pl.close()
