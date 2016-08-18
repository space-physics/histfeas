#!/usr/bin/env python
"""
-m realvid png   # plot simultaneous HST 2 cam + DASC all sky
-m fwd optim png h5 # compute inversion based on optical intensity, save plots as png and result as hdf5 for quick replot

--load -m fwd optim png # load and replot inversion with say different axes limits (far faster than needlessly recomputing)
"""
import logging
logging.basicConfig(level=logging.WARNING)
from dateutil.parser import parse
from sys import argv
#
#from pythonutils.ulimit_nofile import raise_nofile
#raise_nofile(4096) # ulimit -n 1024 will crash with OSError. This is a temporary setting extinguishing with Python session.


def hist_figure(xlsreg,makecomp):
    #imported here to allow matplotlib.use
    from histfeas.main_hist import doSim

    doSim(ParamFN=xlsreg,
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


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('ini',help='ini file to load')
    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-m','--makeplot',help='plots to make (realvid png      fwd optim png h5)',default=[],nargs='+')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',type=int,nargs='+')
    p.add_argument('-o','--outdir',help='output directory',default='out/2013-04-14')
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

    outdir = p.outdir
    timeInds=p.frames

    x1d = None
    vlim = {'p':[None,None,90,400,5e4,8e5,5e4,8e5], 'j':[10,250, 10,250],
            'b':[0,2000]}

    treq = [parse(t) for t in p.treq] if p.treq else None

    overrides = {'treq':treq}

    if not p.load:
        print('writing {} to {}'.format(p.makeplot,outdir))
        hist_figure(p.ini,p.makeplot)
        p.makeplot = [] #don't redo plots just made
#%% load results and plot
    from histfeas.loadAnalyze import readresults,findxlsh5 #import here to allow matplotlib.use
    h5list,inifn = findxlsh5(outdir)
    readresults(h5list,inifn,vlim,x1d,overrides,p.makeplot,p.verbose)

    show()
