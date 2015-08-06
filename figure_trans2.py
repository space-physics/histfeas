#!/usr/bin/env python3
"""
Figure translating generated by HiST program
Michael Hirsch
"""
from sys import argv
from matplotlib.pyplot import show
#
from histfeas.main_hist import doSim
from histfeas.loadAnalyze import readresults,findxlsh5

def hist_figure():
    Phi0,Phifit =doSim(ParamFN=regXLS,
                  makeplot=['fwd','optim','png','h5'],
                  timeInds=timeInds,
                  overrides = overrides, #{'minev': minev,'filter':filt, 'fwdguess':fwdguess,
				                    #'fitm':fitm,'cam':cam,'camx':acx,'ell':ell,'Jfwd':influx},
                  progms = outdir,
                  x1d=x1d,
                  vlim = vlim,
                  animtime=None,
                  cmd = ' '.join(argv),
                  verbose=0
                  )

    return Phi0,Phifit


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='translating figure plotter')
    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-m','--makeplot',help='plots to make',default=[],nargs='+')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p = p.parse_args()

    regXLS='in/2cam_trans.xlsx'
    timeInds=[4,23,47]
    outdir='out/rev2_trans2'
    x1d = [-0.55,1.35,3.75]
    vlim = {'p':[-1.5,4.5,90,300,5e7,8e8,5e7,2e9], 'j':[1e3,1.1e5, 1e3,8e5], 'b':[0,3e3]}
    overrides = {'ell':False}

    if not p.load:
        print('running Hist program -- will write png and h5 to ' + outdir)
        Phi0,Phifit=hist_figure()


    h5list,xlsfn = findxlsh5(outdir)
    readresults(h5list,xlsfn,vlim,x1d,overrides,p.makeplot,p.verbose)

    if 'show' in p.makeplot:
        show()