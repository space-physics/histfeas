#!/usr/bin/env python3
"""
Registration case for HiST program
Michael Hirsch
"""
from numpy import isclose
from numpy.testing import assert_allclose
from sys import argv
import h5py
from tempfile import gettempdir
#from os import devnull
#
from histfeas.main_hist import doSim



def hist_registration(regh5,regXLS):
    Phi0,Phifit =doSim(ParamFN=regXLS,
                  makeplot=['fwd','optim'],
                  timeInds=None,
                  overrides = None, #{'minev': minev,'filter':filt, 'fwdguess':fwdguess, 'fitm':fitm,'cam':cam,'camx':acx,'ell':ell,'Jfwd':influx},
                  progms = gettempdir(),
                  x1d=[None],
                  vlim = {'p':[None]*6,'j':[None]*2,'b':[None]*2},
                  animtime=None,
                  cmd = ' '.join(argv),
                  verbose=0
                  )

    return Phi0,Phifit

def readCheck(Phi0,Phifit):
    with h5py.File(regh5,'r',libver='latest') as f:
        assert_allclose(f['/phifwd/phi'],Phi0[...,0])
        # noise makes inversion result differ uniquely each run
        assert isclose(f['/phifwd/E0'],Phifit[0]['gE0'],rtol=0.2)
        assert isclose(f['/phifwd/x0'],Phifit[0]['gx0'],rtol=0.2)

        #the str() are needed instead of format() !
        print('E0 estimation error [eV] ' +str(f['/phifwd/E0']-Phifit[0]['gE0']))
        print('x0 estimation error [km] ' +str(f['/phifwd/x0']-Phifit[0]['gx0']))


def writeout(regh5):
    with h5py.File(regh5,'a',libver='latest') as f:
        f['/phifwd/E0'] = 7500.
        f['/phifwd/x0'] = 1.

if __name__ == '__main__':
    regh5='test/registration.h5';     regXLS='test/registration.xlsx'

    Phi0,Phifit=hist_registration(regh5,regXLS)
    readCheck(Phi0,Phifit)
