import os

import numpy as np

import pylab as pl
import matplotlib.cm as mlcm
import matplotlib.colors as mlcol
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

from bic import atoms_jumps

def plot_breaks(brk, orientation='v', ls='-'):
    if orientation == 'v':
        axline = pl.axvline
    elif orientation == 'h':
        axline = pl.axhline
    else:
        raise ValueError()

    for b in brk:
        axline(b-0.5, lw=.5, ls=ls, c='k')

def add_chr_ticks(ax, breaks, L):
    chr_ticks = np.array([0] + breaks.tolist() + [L])
    chr_ticks = chr_ticks[:-1] + (chr_ticks[1:] - chr_ticks[:-1])/2.0

    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(['chr%d' % c for c in range(1, len(chr_ticks)+1)],
                       #rotation=45
                       )

def plot_reconstruction(Y, Yest, J, breaks, title):
    fig = pl.figure()
    pl.suptitle('Data reconstruction (%s)' % title, weight='bold', size=8)

    L, S = Y.shape

    # Color boundaries and cmap
    vmin = -2.5 * Y.std()
    vmax = 2.5 * Y.std()
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    cmap = mlcm.RdBu_r

    fit_error = np.sum((Y - Yest)**2.) / (S*L)

    IMAGES = {'$Y$': Y.T, '$\\hat{Y}$ [$Err=%2.3e$]' % fit_error: Yest.T}
    for i, k in enumerate(IMAGES):
        pl.subplot(2,1,i+1)
        im_ref = pl.imshow(IMAGES[k], aspect='auto',
                           interpolation='none', norm=norm, cmap=cmap)

        ax = pl.gca()
        plot_breaks(breaks)
        add_chr_ticks(ax, breaks, L)
        pl.title(k, size=8)

        #plot_breaks(np.arange(1, S, 1), orientation='h', ls='--')
        pl.yticks(range(S), ['S%d' % t for t in range(1, S+1)])
        #pl.yticks(range(19, S, 20), [str(x+1) for x in range(19, S, 20)], size=10)

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(6)
            tick.tick1On = False
            tick.tick2On = False
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(5)
            tick.tick1On = False
            tick.tick2On = False

    pl.tight_layout(pad=2.0)

    # Make an axis for the colorbar on the right side
    pl.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cb = fig.colorbar(im_ref, cax=cax, extend='both')
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(6)

def plot_representation(Theta, breaks, atoms_order):
    fig = pl.figure()
    pl.suptitle(r'Representation coefficients ($\Theta$)',
                weight='bold', size=8)

    J, S = Theta.shape

    # Color boundaries and cmap
    tmin = Theta.min()
    vmin = 0.0 if tmin >= 0.0 else -np.ceil(np.abs(tmin))
    vmax = np.ceil(Theta.max())
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    if vmin == 0.0:
        cmap = mlcm.Reds
    else:
        cmap = mlcm.RdBu_r

    im_ref = pl.imshow(Theta[atoms_order], aspect='auto',
                       interpolation='none', norm=norm, cmap=cmap)
    plot_breaks(np.arange(1, S, 1), orientation='v', ls='--')
    plot_breaks(np.arange(1, J, 1), orientation='h')

    pl.xticks(range(S), ['S%d' % t for t in range(1, S+1)], size=5)
    #pl.xticks(range(19, S, 20), [str(x+1) for x in range(19, S, 20)], size=10)
    pl.yticks(range(J), ['A%d' % (t+1) for t in atoms_order], size=10)

    cb = pl.colorbar(extend='both', orientation='horizontal')
    #for t in cb.ax.get_xticklabels():
    #    t.set_fontsize(5)

def plot_atoms(B, breaks, atoms_order, epsj=1e-3, gth=None, lth=None):
    C = 2
    R = len(atoms_order)/C
    if R*C < len(atoms_order):
        R += 1

    fig = pl.figure(figsize=(10, R))

    jumps = atoms_jumps(B, epsj)
    pl.suptitle('Dictionary Atoms ($B$, jumps=$%d$)' % jumps,
                weight='bold', size=8)

    L, J = B.shape

    # Color boundaries and cmap
    vmax = np.abs(B).max()
    vmin = -vmax
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    cmap = mlcm.RdBu_r

    #division = np.sqrt(len(atoms_order))
    #R, C = int(np.ceil(division)), int(np.floor(division))

    for p, d in enumerate(B[:,atoms_order].T):
        pl.subplot(R, C, p+1)

        pl.plot(d, 'k-', lw=1, alpha=0.5)
        im_ref = pl.scatter(range(len(d)), d, c=d, s=40, marker='.',
                            edgecolors='none', cmap=cmap, norm=norm)

        vbound = 1.4 #np.abs(d).max()*1.1
        pl.axis([0, len(d), -vbound, vbound])

        ax = pl.gca()

        plot_breaks(breaks)
        add_chr_ticks(ax, breaks, L)
        ax.yaxis.set_label_position('right')
        ax.set_ylabel('Atom #%d' % (atoms_order[p]+1),
                      rotation=-90, size=6, weight='bold')
        
        if gth:
            pl.axhline(gth, c='r')
        if lth:
            pl.axhline(lth, c='c')

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(6)
            tick.tick1On = False
            tick.tick2On = False
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(6)

    # Spacing
    fig.tight_layout()
    pl.subplots_adjust(top=0.85)
    
    # Make an axis for the colorbar on the right side
    pl.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cb = fig.colorbar(im_ref, cax=cax, extend='both')
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(6)

# Multipage pdf with plots --------------------------------------------
def pdf_plots(filename, Y, Theta, B, J, lambda_, mu, tau, w, epsj):
    pp = PdfPages(filename)

    atoms_order = np.argsort(Theta.sum(axis=1))[::-1]
    breaks = (np.where(w==0)[0]) if not w is None else []

    # Parameters printed into the data reconstruction title
    parstr = ('$J={:2d}$, $\\lambda={:10.5f}$, '
              '$\\mu={:10.5f}$, '
              '$\\tau={:10.5f}$').format(J, lambda_, mu, tau)

    plot_reconstruction(Y, np.dot(B, Theta),
                        J, breaks, title=parstr)
    pp.savefig()

    plot_representation(Theta, breaks, atoms_order)
    pp.savefig()

    plot_atoms(B, breaks, atoms_order, epsj)
    pp.savefig()

    pp.close()
    del pp

def png_plots(filename, Y, Theta, B, J, lambda_, mu, tau, w, epsj):
    atoms_order = np.argsort(Theta.sum(axis=1))[::-1]
    breaks = (np.where(w==0)[0]) if not w is None else []

    # Parameters printed into the data reconstruction title
    parstr = ('$J={:2d}$, $\\lambda={:10.5f}$, '
              '$\\mu={:10.5f}$, '
              '$\\tau={:10.5f}$').format(J, lambda_, mu, tau)

    pl.figure()
    plot_reconstruction(Y, np.dot(B, Theta),
                        J, breaks, title=parstr)
    pl.savefig('%s_reconstruction.png' % filename)
    pl.clf()
    pl.close()

    pl.figure()
    plot_representation(Theta, breaks, atoms_order)
    pl.savefig('%s_representation.png' % filename)
    pl.clf()
    pl.close()

    pl.figure()
    plot_atoms(B, breaks, atoms_order, epsj)
    pl.savefig('%s_atoms.png' % filename)
    pl.clf()
    pl.close()


def plot_bics(filename, J_range, lambda_range, mu_range, tau_range, BICS):
    fig = pl.figure()
    
    mu_range = sorted(mu_range)
    tau_range = sorted(tau_range)

    # BICS collected in (mu, lambda, tau) order
    m, l, t = len(mu_range), len(lambda_range), len(tau_range)
    X, Y = np.meshgrid(mu_range, tau_range, indexing='ij')

    cm = pl.get_cmap('gist_rainbow')
    colors = [cm(x) for x in np.linspace(0, 1, len(J_range)*l)]

    idx = 0
    ax = fig.add_subplot(111, projection='3d')
    for j, J in enumerate(J_range):
        Z = np.array(BICS[J]).reshape((m, l, t))
        
        for i, lam in enumerate(sorted(lambda_range)):
            title = r'$J=%d$, $\lambda=%2.3e$' % (J, lam)
            im_ref = ax.plot_wireframe(X, Y, Z[:,i,:], rstride=10, cstride=10,
                                       linewidth=1, antialiased=False,
                                       label=title, color=colors[idx])            
            idx+=1

            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(6)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(6)
            for tick in ax.zaxis.get_major_ticks():
                tick.label1.set_fontsize(6)

            ax.set_xlabel(r'$\mu$')
            ax.set_ylabel(r'$\tau$')
            ax.set_zlabel('$BIC$')

    pl.legend(loc='upper left', ncol=2, fancybox=True, prop={'size':8})
    pl.tight_layout()
    pl.savefig(filename)
