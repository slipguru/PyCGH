# Author: Salvatore Masecchia <salvatore.masecchia@disi.unige.it>
# License: New BSD
""" Optimization main functions """

import os
import sys
import logging
import datetime as dt

import numpy as np

from bic import pseudo_BIC_search as BIC_search
   
### Imports current parameters ------------------------------------------------
from config import *
    
if __name__ =='__main__':
    ### Data ------------------------------------------------------------------
    data = np.load(sys.argv[1])
    Y = data['Y_0.1']
    #Y = data['Y_1.0']
    chrbrks = data['chrbreaks']
    if params['tvw'] == 'auto':
        params['tvw'] = chrbrks
    
    ### Output ----------------------------------------------------------------
    OUTPUT_DIR = '%s/Results' % OUTPUT
    
    ### Logging and saving ----------------------------------------------------
    logger = logging.getLogger('cghDL')
    logger.setLevel(logging.DEBUG)
    logfilename = 'TIME_%s_%s.log' % (ALGORITHM, dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    lfile = logging.FileHandler(os.path.join(OUTPUT_DIR, logfilename))
    lfile.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    lfile.setFormatter(formatter)
    logger.addHandler(lfile)

    infomsg = '* Experiment with J={:2d}, lambda={:10.5f}, mu={:10.5f}, tau={:10.5f}'
    convmsg = ('  - time: {:<10.3f}s\n')

    # CALL-BACK FUNCTION: called after each CGHDL fit #########################
    def callback(time, J, lambda_, mu, tau):
        msg = infomsg.format(J, lambda_, mu, tau)
        logger.info(msg)
        msg = convmsg.format(time)
        logger.info('\n'+msg)
   ###########################################################################
    
    params['callback'] = callback
    result = BIC_search(Y, **params)
        
    