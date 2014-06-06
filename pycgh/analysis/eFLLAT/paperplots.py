import os
import sys

import numpy as np

import pylab as pl
import matplotlib.cm as mlcm
import matplotlib.colors as mlcol
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

FONT_SIZE = 3
DPI = 600

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
    ax.set_xticklabels(['chr%d' % c for c in range(1, len(chr_ticks)+1)], size=FONT_SIZE)

def plot_atoms(B, breaks, atoms_order, positions, grid):
    L, J = B.shape
    
    # Color boundaries and cmap
    vmax = 1.0
    vmin = -1.0
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    cmap = mlcm.RdBu_r

    for p, (d, pos) in enumerate(zip(B[:,atoms_order].T, positions)):
        pl.subplot2grid(grid, pos)

        pl.plot(d, 'k-', lw=1, alpha=0.5)
        im_ref = pl.scatter(range(len(d)), d, c=d, s=40, marker='.',
                            edgecolors='none', cmap=cmap, norm=norm)

        vbound = 1.4 #np.abs(d).max()*1.1
        pl.axis([0, len(d), -vbound, vbound])

        ax = pl.gca()

        plot_breaks(breaks)
        add_chr_ticks(ax, breaks, L)
        ax.yaxis.set_label_position('right')
        ax.set_ylabel('Atom #%d' % (p+1), rotation=-90, size=FONT_SIZE, weight='bold')
        
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(FONT_SIZE)
            tick.tick1On = False
            tick.tick2On = False
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(FONT_SIZE)
        
    # FAKE!
    for p, (d, pos) in enumerate(zip(B[:,atoms_order].T, [(0, 0), (0, 1)])):
        pl.subplot2grid(grid, pos)

        pl.plot(d, 'k-', lw=1, alpha=0.5)
        im_ref = pl.scatter(range(len(d)), d, c=d, s=40, marker='.',
                            edgecolors='none', cmap=cmap, norm=norm)

        vbound = 1.4 #np.abs(d).max()*1.1
        pl.axis([0, len(d), -vbound, vbound])

        ax = pl.gca()

        plot_breaks(breaks)
        add_chr_ticks(ax, breaks, L)
        ax.yaxis.set_label_position('right')
        ax.set_ylabel('Atom FAKE #%d' % (p+1), rotation=-90, size=FONT_SIZE, weight='bold')
        
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(FONT_SIZE)
            tick.tick1On = False
            tick.tick2On = False
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(FONT_SIZE)   
        
def plot_representation(Theta, breaks, atoms_order, grid):
    J, S = Theta.shape
       
    pl.subplot2grid(grid, (grid[0]-6, 0), rowspan=2, colspan=2)

    # Color boundaries and cmap
    vmax = 1.0
    vmin = -1.0
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    cmap = mlcm.RdBu_r

    im_ref = pl.imshow(Theta[atoms_order], aspect='auto',
                       interpolation='none', norm=norm, cmap=cmap)
    plot_breaks(np.arange(1, S, 1), orientation='v', ls='--')
    plot_breaks(np.arange(1, J, 1), orientation='h')

    pl.xticks(range(S), ['S%d' % t for t in range(1, S+1)])
    pl.yticks(range(J), ['A#%d' % t for t in range(1, J+1)], weight='bold')
    
    ax = pl.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE)
        tick.tick1On = False
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE)
        tick.tick1On = False
        tick.tick2On = False

def plot_reconstruction(Y, fit_error, J, breaks, grid):
    L, S = Y.shape

    pl.subplot2grid(grid, (grid[0]-4, 0), rowspan=4, colspan=2)

    # Color boundaries and cmap
    vmin = -1.0
    vmax = 1.0
    norm = mlcol.normalize(vmax=vmax, vmin=vmin)
    cmap = mlcm.RdBu_r

    if fit_error:
        error_tit = r'\||\mathbf{\bar{Y}} - \mathbf{\hat{Y}}\||^2 / (S \cdot L) = %2.3e' % fit_error
        title = r'$\mathbf{\hat{Y}} \quad\quad [ %s ] $' % error_tit    
    else:
        title = r'$\mathbf{Y}$'
    im_ref = pl.imshow(Y.T, aspect='auto', interpolation='none',
                       norm=norm, cmap=cmap)

    ax = pl.gca()
    plot_breaks(breaks)
    add_chr_ticks(ax, breaks, L)
    pl.title(title, size=FONT_SIZE+2)

    pl.yticks(range(S), ['S%d' % t for t in range(1, S+1)])

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE)
        tick.tick1On = False
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE)
        tick.tick1On = False
        tick.tick2On = False
    return ax

def plot_dictionary(B, breaks, atoms_order):
        
    # Dictionary plotting
    pl.figure(figsize=(4.25, 7))
    ax = pl.subplot(211)
    
    for idx, atom in enumerate(B[:,atoms_order].T):
        ax.semilogy(np.abs(atom), '-', lw=1, label='Atom #%d'%(idx+1))
    plot_breaks(breaks)
    add_chr_ticks(ax, breaks, B.shape[0])
    pl.axis([0, B.shape[0], 1e-18, 1e1])
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE+2)
        tick.tick1On = False
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE+2)
        tick.tick1On = False
        tick.tick2On = False
       
    ax = pl.subplot(212)
       
    for idx, atom in enumerate(B[:,atoms_order].T):
        ax.semilogy(np.abs(np.diff(atom)), '-', lw=1, label='Atom #%d'%(idx+1))
    plot_breaks(breaks)
    add_chr_ticks(ax, breaks, B.shape[0])
    pl.axis([0, B.shape[0], 1e-18, 1e1])
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE+2)
        tick.tick1On = False
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(FONT_SIZE+2)
        tick.tick1On = False
        tick.tick2On = False
    
    pl.tight_layout()
    
    # Shink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=4, prop={'size':FONT_SIZE+4})
      

def summary_plot(filename, Y, Theta, B, J, w, positions, atoms_order=None, fit_error=None):
    atoms_order = range(J) if atoms_order is None else atoms_order
    breaks = (np.where(w==0)[0]) if not w is None else []

    def cm2inch(value):
        return value/2.54

    #if len(positions) == 5: # original data
    #            # larghezza,  altezza
    #    #size = (cm2inch(15), cm2inch(20))
    #    size = (4.25, 5)
    #    grid = (9, 2)
    #elif len(positions) == 7: #results
    #            # larghezza,  altezza
    #    #size = (cm2inch(15), cm2inch(22))
    size = (4.25, 7)
    grid = (10, 2)
    #else:
    #    raise Exception('wrong positions!')

    fig = pl.figure(figsize=size)
    pl.subplot2grid(grid,(0, 0))    
    plot_atoms(B, breaks, atoms_order, positions, grid)
    plot_representation(Theta, breaks, atoms_order, grid)  
    ax = plot_reconstruction(Y, fit_error, J, breaks, grid)  # no fit error
    
    cb = pl.colorbar(extend='both', orientation='horizontal', use_gridspec=True)
    for t in cb.ax.get_xticklabels():
        t.set_fontsize(FONT_SIZE)
    
    pl.tight_layout(h_pad=0.0)
    pl.savefig('%s_atoms_repres_reconstr.png' % filename, dpi=DPI)
    pl.clf()
    pl.close()
    
    plot_dictionary(B, breaks, atoms_order)
    pl.savefig('%s_dictionary.png' % filename, dpi=DPI)
    pl.clf()
    pl.close()

if __name__ == '__main__':
    data = np.load(sys.argv[1])
    w = data['chrbreaks']
    
    if sys.argv[2] == 'NONE':
        L, S, J = 500, 25, 5
        Y0p1 = data['Y_0.1']
        Y = data['YP']          # con atomi aggiunti
        YT = data['Y']          # senza atomi aggiunti
        Theta = data['Theta']
        B = data['B']
          
        filename = os.path.split(os.path.splitext(sys.argv[1])[0])[-1]  
        summary_plot(os.path.join(sys.argv[3], filename), Y0p1, Theta, B, J, w,
                     positions=[(1, 0), (1, 1), (2, 0), (2, 1), (3, 0)])
    else:
        sol = np.load(sys.argv[2])
        B = sol['B']
        Theta = sol['Theta']
        Y = data['YP']
        Yest = np.dot(B, Theta)

        Y[:,-5:] = 0.0                
        
        try:
            J = sol['J']
            lambda_ = sol['lambda']
            mu = sol['mu']
            tau = sol['tau']
        except KeyError:
            J = 7
            lambda_ = 1e-4
            mu = 1.0
            tau = 1e-3
            print '*******Manual parameters!!******'
        
        S, L = Y.shape
        fit_error = np.sum((Y - Yest)**2.) / (S*L)
    
        filename = 'J_%02d-lambda_%02.3e-mu_%02.3e-tau_%02.3e' % (J, lambda_,
                                                                  mu, tau)
        summary_plot(os.path.join(sys.argv[3], filename),
                     Yest, Theta, B, J, w,
                     positions=[(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0)],
                     atoms_order = np.array([6,3,5,1,7,4,2]) - 1,
                     fit_error=fit_error)
            

    
    