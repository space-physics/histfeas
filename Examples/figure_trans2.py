#!/usr/bin/env python
"""
Figure translating generated by HiST program
"""

#
from histfeas.main_hist import doSim
from histfeas.loadAnalyze import readresults,findxlsh5

def hist_figure(xlsreg):
    Phi0,Phifit =doSim(ParamFN=xlsreg,
                  makeplot=['fwd','optim','png','h5'],
                  timeInds=timeInds,
                  overrides = overrides, #{'minev': minev,'filter':filt, 'fwdguess':fwdguess,
				                    #'fitm':fitm,'cam':cam,'camx':acx,'ell':ell,'Jfwd':influx},
                  odir = outdir,
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
    p.add_argument('--ell',help='compute projection matrix',action='store_true')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',type=int,default=(1,3))
    p.add_argument('-o','--outdir',help='output directory',default='out/rev2_trans2')
    p = p.parse_args()

    xlsreg='in/2cam_trans.xlsx'
    outdir = Path(p.outdir).expanduser()
    timeInds=p.frames

    x1d = [1.55,3.75]
    vlim = {'p':[5e7,8e8],  'p1d':(5e7,2e9),
            'j':[1e3,1.1e5],'j1d':(1e3,8e5),
            'b':[0,1.5e3],
            'x':(-1.5,4.5), 'z':(90,300)}
    overrides = {'ell':p.ell}

    if not p.load:
        print('running HiSTfeas program -- will write png and h5 to {}'.format(outdir))
        Phi0,Phifit=hist_figure(xlsreg)

    h5list,xlsfn = findxlsh5(outdir)
    readresults(h5list,xlsfn,vlim,x1d,overrides,p.makeplot,p.verbose)

    if 'show' in p.makeplot:
        show()