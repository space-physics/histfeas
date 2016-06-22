#!/usr/bin/env python
"""
Registration case for HiST program
Michael Hirsch

To generate registration.h5 (only when making a new type of sim or real) type
python
"""
from histfeas import Path
from os import chdir
from numpy.testing import assert_allclose
from sys import argv
import h5py
from tempfile import gettempdir
#
import matplotlib
matplotlib.use('Agg') #fixes Travis NO DISPLAY bug
#
from histfeas.main_hist import doSim

rootdir = Path(__file__).parents[1]

def hist_registration(regh5,regXLS,odir):
    """
    This creates output hdf5 files, the typical output of the program for offline analysis
    """

    doSim(ParamFN=regXLS,
          makeplot=['fwd','optim','h5'],
          timeInds=None,

          overrides = {'rootdir':rootdir},#,'ell':True}, #{'minev': minev,'filter':filt, 'fwdguess':fwdguess, 'fitm':fitm,'cam':cam,'camx':acx,'ell':ell,'Jfwd':influx},
          odir = odir,
          x1d=None,
          vlim = {'p':[None]*6,'j':[None]*2,'b':[None]*2},
          animtime=None,
          cmd = ' '.join(argv),
          verbose=0
          )


def readCheck(Phi0,Phifit):
    with h5py.File(str(regh5),'r',libver='latest') as f:
        assert_allclose(f['/phifwd/phi'],Phi0[...,0])
        # noise makes inversion result differ uniquely each run
        xerrpct=(f['/phifwd/x0'] - Phifit[0]['gx0']) / f['/phifwd/x0'] * 100
        xmsg = 'x0 estimation error [km] {:.1f} %'.format(xerrpct[0])
        assert abs(xerrpct[0]) < 30,'B_\perp location error out of tolerance'

        Eerrpct = (f['/phifwd/E0'] - Phifit[0]['gE0']) / f['/phifwd/E0'] * 100
        Emsg = 'E0 estimation error [eV] {:.1f} %'.format(Eerrpct[0])
        assert abs(Eerrpct[0]) < 30,'E_0 error out of tolerance'

        print(xmsg)
        print(Emsg)

def writeout(regh5):
    with h5py.File(str(regh5),'a',libver='latest') as f:
        f['/phifwd/E0'] = 7500.
        f['/phifwd/x0'] = 1.

if __name__ == '__main__':
    from histfeas.loadAnalyze import readresults,findxlsh5 #here for matplotlib import

    chdir(str(rootdir)) # in case running from other than project root directory.
#%% simulation only
    tdir  = Path(__file__).parent
    regh5 = tdir / 'registration.h5'
    regini= tdir / 'registration.ini'

    odir=gettempdir()
#%% do inversion
    hist_registration(regh5,regini,odir)
#%% find result HDF5
    h5list,_ = findxlsh5(odir)
#%% load result
    Phi0,Phifit = readresults(h5list,regini)
#%% check vs known result
    readCheck(Phi0,Phifit)
    print('OK:  simulation registration case')
#%% real data
