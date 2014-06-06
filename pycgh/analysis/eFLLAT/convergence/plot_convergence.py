import sys
import os
import pylab as pl
import numpy as np

result = np.load(sys.argv[1])
path = os.path.dirname(sys.argv[1])
filename = os.path.join(path, 'Result')

PSI_GAPS = result['PSI_GAPS']
PHI_GAPS = result['PHI_GAPS']

fig = pl.figure(figsize=(10, 10))
pl.suptitle(r'PSI_GAPS', weight='bold', size=8)
for p, gaps in enumerate(PSI_GAPS):
    pl.subplot(5, 2, p+1)
    pl.semilogy(gaps[0], '-b')
    pl.semilogy(gaps[1], '-r')
    pl.axhline(np.min(gaps), c='k', lw=1, label='%.3e' % np.min(gaps))
    pl.legend(loc='upper right')
fig.tight_layout()
pl.subplots_adjust(top=0.95)
pl.savefig('%s_LogY_GAPS_PSI.png' % filename)

fig = pl.figure(figsize=(10, 10))
pl.suptitle(r'PHI_GAPS', weight='bold', size=8)
for p, gaps in enumerate(PHI_GAPS):
    pl.subplot(5, 2, p+1)
    pl.semilogy(gaps[0], '-b')
    pl.semilogy(gaps[1], '-r')
    pl.axhline(np.min(gaps), c='k', lw=1, label='%.3e' % np.min(gaps))
    pl.legend(loc='upper right')
fig.tight_layout()
pl.subplots_adjust(top=0.95)
pl.legend()
pl.savefig('%s_LogY_GAPS_PHI.png' % filename)
    
#pl.show() 
