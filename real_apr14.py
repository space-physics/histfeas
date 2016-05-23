#!/usr/bin/env python3
"""
-c realvid png   # plot simultaneous HST 2 cam + DASC all sky
-c fwd optim png h5 # compute inversion based on optical intensity, save plots as png and result as hdf5 for quick replot

"""
import logging
logging.basicConfig(level=logging.WARNING)
from tempfile import gettempdir
from pathlib import Path
from dateutil.parser import parse
from sys import argv

def hist_figure(xlsreg,makecomp):
    #imported here to allow matplotlib.use
    from histfeas.main_hist import doSim

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
    p.add_argument('-m','--makeplot',help='plots to make (realvid png      fwd optim png h5)',default=[],nargs='+')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',type=int,nargs='+')
    p.add_argument('-o','--outdir',help='output directory',default=Path(gettempdir()) / 'out/apr14T085430')
    p.add_argument('-t','--treq',help='specific times requested',nargs='+')
    p = p.parse_args()

    if not 'show' in p.makeplot:
        import matplotlib
        matplotlib.use('Agg')
        print(matplotlib.get_backend())
    #
    from matplotlib.pyplot import show

    import seaborn as sns
    sns.color_palette("cubehelix")
    sns.set(context='paper', style='whitegrid',font_scale=1.5,
        rc={'image.cmap': 'cubehelix_r'})


    xlsreg='in/apr14T085430.ini'
    outdir = Path(p.outdir)
    timeInds=p.frames

    x1d = None
    vlim = {'p':[-3,7,90,400,5e4,8e5,5e4,8e5], 'j':[10,250, 10,250],
            'b':[0,2000]}

    treq = [parse(t) for t in p.treq] if p.treq else None

    overrides = {'treq':treq}

    if not p.load:
        print('running HiSTfeas program writing {} to {}'.format(p.makeplot,outdir))
        Phi0,Phifit = hist_figure(xlsreg,p.makeplot)

    from histfeas.loadAnalyze import readresults,findxlsh5 #import here to allow matplotlib.use
    h5list,xlsfn = findxlsh5(outdir)
    readresults(h5list,xlsfn,vlim,x1d,overrides,p.makeplot,p.verbose)

    if 'show' in p.makeplot:
        show()
