#!/usr/bin/env python3
"""
Registration case for HiST program
Michael Hirsch
"""
from numpy import isclose
from numpy.testing import assert_allclose
from sys import argv
import h5py
#from os import devnull
#
from main_hist import doSim

regh5='test/registration.h5'

def hist_registration():
    Phi0,Phifit =doSim(ParamFN='test/registration.xlsx',
                  makeplot=['fwd','optim'],
                  timeInds=None,
                  overrides = None, #{'minev': minev,'filter':filt, 'fwdguess':fwdguess, 'fitm':fitm,'cam':cam,'camx':acx,'ell':ell,'Jfwd':influx},
                  progms = None,
                  x1d=[None],
                  vlim = {'p':[None]*6,'j':[None]*2,'b':[None]*2},
                  animtime=None,
                  cmd = ' '.join(argv),
                  verbose=0
                  )

    return Phi0,Phifit

def readCheck(Phi0,Phifit):
    with h5py.File(regh5,'r',libver='latest') as f:
        assert_allclose(f['/phifwd/phi'],Phi0)
        # noise makes inversion result differ uniquely each run
        assert isclose(f['/phifwd/E0'],Phifit[0]['gE0'],rtol=0.15)
        assert isclose(f['/phifwd/x0'],Phifit[0]['gx0'],rtol=0.15)

        print('E0 estimation error is {} eV  {} km'.format(
                                    f['/phifwd/E0']-Phifit[0]['gE0'],
                                    f['/phifwd/x0']-Phifit[0]['gx0']
                                    ))

def writeout(Phi0,Phifit):
    with h5py.File(regh5,'w',libver='latest') as f:
        f['/phifwd/phi'] = Phi0
        f['/phifwd/E0'] = 6687.
        f['/phifwd/x0'] = 1.
        f['/phifwd/Ek'] = Phifit[0]['EK']

if __name__ == '__main__':
    Phi0,Phifit=hist_registration()

    readCheck(Phi0,Phifit)
