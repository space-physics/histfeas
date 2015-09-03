#!/usr/bin/env python3
"""
Example command line interfacer for HIST feasibility
Michael Hirsch
"""
from signal import signal,SIGINT #for Ctrl C
from os.path import expanduser
from numpy import arange
import matplotlib as mpl
from matplotlib.pyplot import show
from sys import argv
#
from histfeas.main_hist import doSim

def signal_handler(signal, frame):
    exit('\n *** Aborting program as per user pressed Ctrl+C ! \n')

if __name__ == '__main__':
    from argparse import ArgumentParser

    signal(SIGINT, signal_handler) #allows friendly ctrl-c interrupt

    p = ArgumentParser(description='analyzes HST data and makes simulations')
    p.add_argument('infile',help='.xls filename with simulation parameters')
    p.add_argument('outdir',help='directory for output')
    p.add_argument('-v','--verbose',help='set debugging verbosity e.g. -v -vv -vvv',action='count',default=0)
    p.add_argument('-m','--makeplot',help='list plots you want made',nargs='+',default=['']) #None gave type errors in generators & list comp.
    p.add_argument('-f','--frames',help='START STOP STEP of time indices',nargs=3,type=int)
    p.add_argument('--profile',help='Profile performance of program (development/debug only)',action='store_true')
    p.add_argument('-c','--cam',help='zero-indexed cameras to use (overrides XLS)',nargs='+',default=[None],type=int)
    p.add_argument('--cx',help='override cam positions (must specify all)',nargs='+',type=float)
    p.add_argument('--influx',help='flux .h5 file to use (overrides XLS)')
    p.add_argument('--x1d',help='required location [km] of x for 1-D flux plots',nargs='+',default=[None],type=float)
#    p.add_argument('--saveall',help='saves "all" variables to disk, useful for registration case tracking',action='store_true')
    p.add_argument('-a','--anim',help='animate plots (crudely)',type=float)
    p.add_argument('--filter',help='optical filter choices: bg3   none')
    p.add_argument('--minev',help='minimum beam energy to include in fwd model',type=float)
    p.add_argument('-g','--fwdguess',help='feed minimizer fwd answer. true | randn stddev |',nargs='+',default=[None])
    p.add_argument('--fitm',help='override fit (minimization) method')
    p.add_argument('--vlim',help='xmin xmax zmin zmax pmin pmax   limits for VER plots',type=float,nargs=6,default=[None]*6)
    p.add_argument('--jlim',help='MIN MAX flux limits for diff num flux plots',type=float,nargs=2,default=(None,)*2)
    p.add_argument('--blim',help='MIN MAX flux limits for brightness plots',type=float,nargs=2,default=(None,)*2)
    p.add_argument('-L','--ell',help='force recomputation of sparse matrix L',action='store_true')
    ar = p.parse_args()

    ParamFN = expanduser(ar.infile)
#%% plot setup
    makeplot = ar.makeplot
    print(('matplotlib backend / version: ' + mpl.get_backend() +'  ' + mpl.__version__  ))
#%%
    vlim = {'p':ar.vlim,'j':ar.jlim,'b':ar.blim}
#%%
    timeInds = ar.frames
    doProfile = ar.profile
#%% overrides
    overrides = {'minev': ar.minev,'filter':ar.filter,
                 'fwdguess':ar.fwdguess, 'fitm':ar.fitm,'cam':ar.cam,
                 'camx':ar.cx,'ell':ar.ell,'Jfwd':ar.influx}

    if timeInds is not None:
        timeInds = arange(timeInds[0],timeInds[1],timeInds[2]) #NOT range!!

    if doProfile: #devel debug
        #from subprocess import Popen, PIPE
        import cProfile,pstats
        proffn = 'hstprofstats.pstats'
        print('saving profile results to ' + proffn)
        cProfile.run('doSim(ParamFN,makeplot,timeInds,'
                 'overrides,ar.outdir,ar.x1d,vlim,ar.anim,' '.join(argv),ar.verbose)',proffn)
        pstats.Stats(proffn).sort_stats('time','cumulative').print_stats(50)
        #binpath = expanduser('~/profile/')
        #sysCall = [binpath + 'gprof2dot.py','-f','pstats',profFN,'|','dot','-Tpng','-o','output.png']
        #print(sysCall)
        #po = Popen(sysCall, stdout=PIPE, cwd=binpath, shell=False)
        #so,serr = po.communicate(timeout=1) #timeout is incompatible with Python 2.
        #print(so.decode('utf8'))

    else: #normal
        doSim(ParamFN,makeplot,timeInds,overrides,ar.outdir,ar.x1d, vlim,ar.anim,' '.join(argv), ar.verbose)

    if 'show' in makeplot:
        show()
