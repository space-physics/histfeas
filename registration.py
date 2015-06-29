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

def hist_registration():
    Phi0,Phifit =doSim(ParamFN='in/registration.xlsx',
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

    with h5py.File('test/registration.h5','r',libver='latest') as f:
        assert_allclose(f['/Phi/0'],Phi0)
        # noise makes inversion result differ uniquely each run
        # here I allow 10% error to pass
        assert isclose(f['/Phi/params/E0'],Phifit[0]['gE0'],rtol=0.1)
        assert isclose(f['/Phi/params/x0'],Phifit[0]['gx0'],rtol=0.1)

    return Phi0,Phifit

if __name__ == '__main__':
    Phi0,Phifit=hist_registration()