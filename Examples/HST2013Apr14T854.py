#!/usr/bin/env python
"""
Load real data from HiST only
Michael Hirsch
"""
from histfeas import Path
from dateutil.parser import parse
from pandas import date_range
from sys import argv
from matplotlib.pyplot import show
#
import seaborn as sns
sns.color_palette("cubehelix")
sns.set(context='paper', style='whitegrid',font_scale=1.5,#2,
        rc={'image.cmap': 'cubehelix_r'})
#
from histfeas.main_hist import doSim
from histfeas.loadAnalyze import readresults,findxlsh5
#
import histfeas
root = Path(histfeas.__file__).parents[1]

def hist_figure(xlsreg,makecomp):
    Phi0,Phifit =doSim(ParamFN=xlsreg,
                  makeplot=makecomp,
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
    p = ArgumentParser(description='flaming figure plotter')
    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-c','--compute',help='plot modes to compute the first time  e.g. "realvid png"  or  "optim png h5"',nargs='+',default=['fwd','optim','realvid','png','h5'])
    p.add_argument('-m','--makeplot',help='plots to load and make (must have been first computed via -c)',default=[],nargs='+')
    p.add_argument('--ell',help='compute projection matrix',action='store_true')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',type=int,nargs='+')
    p.add_argument('-o','--outdir',help='output directory',default='~/data/out/realtry')
    p.add_argument('-t','--treq',help='specific times requested',nargs='+')
    p = p.parse_args()

    xlsreg='HST/apr14.xlsx'
    outdir = Path(p.outdir).expanduser()
    timeInds=p.frames
    x1d = None
    vlim = {'p':[-5,10,90,400,5e3,5e5], 'j':[10,100],
            'b':[0,2000]}

    if p.treq:
        treq = [parse(t) for t in p.treq]
    else:
        treq = date_range(start='2013-04-14T8:54:21Z',
                          end=  '2013-04-14T8:54:31Z',
                          freq='1S').to_pydatetime()



    overrides = {'ell':p.ell,
                 'treq':treq,
                 'rootdir':root}

    if not p.load:
        print('running HiSTfeas program -- will write png and h5 to {}'.format(outdir))
        Phi0,Phifit = hist_figure(xlsreg,p.compute)

    h5list,xlsfn = findxlsh5(outdir)
    readresults(h5list,xlsfn,vlim,x1d,overrides,p.makeplot,p.verbose)

    if 'show' in p.makeplot:
        show()
