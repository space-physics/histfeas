#!/usr/bin/env python3
"""
 Michael Hirsch
 May 2014
 GPL v3+

 Procedure:
 0) ingest simulation parameters from XLS file, sanity check, load data
 1) Make an auroral phantom (if not real data, a Simulation, and optionally for a first guess for reconstruction)
 2) create projection matrix
 3) fit differential number flux to reconstructed VER
 4) Test reconstruction & fit by reprojecting to synthetic image brightness

example: (fwd model only)
python3 main_hist.py in/2cam_trans.xlsx /dev/shm/rev_trans2/ -m fwd png --vlim -0.5 3.5 90 350 1e9 1e10 --jlim 1e5 5e5 --blim 0 1e4 -f 0 120 20
"""
from __future__ import division,print_function
from signal import signal,SIGINT #for Ctrl C
from sys import argv
from os.path import expanduser, join
from os import makedirs
from numpy import absolute,zeros,in1d,arange,outer
from numpy.random import normal
from warnings import warn
#
#import logging
#logging.basicConfig(filename='hist.log',filemode='w',level=logging.DEBUG)

def doSim(ParamFN,makeplot,timeInds,overrides,progms,x1d,vlim,animtime, cmd,verbose):
    # local -- these were put here so that matplotlib backend autoselect could happen first
    from matplotlib.pyplot import close
    from pyimagevideo.imageconv import png2multipage
    from sanityCheck import getParams
    from gridaurora.eFluxGen import maxwellian
    from AuroraFwdModel import getSimVER
    from transcararc import getMp,getPhi0,getpx #calls matplotlib
    from observeVolume import getEll,getObs #calls matplotlib
    from FitVER import FitVERopt as FitVER #calls matplotlib
#    from simulFrame import getSimulData
    from plotsnew import goPlot
    from analysehst import analyseres
#%% Step 0) load data
    ap,sim,cam,Fwd = getParams(ParamFN, overrides,makeplot,progms,verbose)
#%% setup loop
    if sim.realdata:
        cam,rawdata,sim = getSimulData(sim,cam,makeplot,progms,verbose)
        sim.nTimeSliceReq = cam['0'].keo.shape[1] #FIXME assumes equal num. of time slices list(cam)[0]
    else: #simulation
        rawdata = None
        if sim.raymap == 'astrometry': #and any('b' in m[:2] for m in makeplot):
            from get1Dcut import get1Dcut #we need cam.angle_deg for plotting
            cam = get1Dcut(cam,makeplot,progms,verbose)
    timeInds = sim.maketind(timeInds)
#%% Step 1) get projection matrix
    Lfwd,Fwd,cam = getEll(sim,cam,Fwd,makeplot,verbose)
#%% preallocation
    PhifitAll = []; drnAll = []; bfitAll=[]
#%% load eigenprofiles from Transcar
    Peig = getMp(sim,Fwd['z'],makeplot,verbose)
#%% synthetic diff. num flux
    if not sim.realdata:
        Phi0all = getPhi0(sim,ap,Fwd['x'],Peig['Ek'], makeplot,verbose)
    if verbose>0: print('timeInds: {}'.format(timeInds))
#%%start looping for each time slice in keogram (just once if simulated)
    for ti in timeInds:
        if sim.realdata:
            Phi0 = None; Pfwd = None
        else: #sim
            """
            we need to integrate in time over the relevant time slices
            the .sum(axis=2) does the integration/smearing in time
            """
            Phi0 = Phi0all[...,ti]
#%% Step 1) Forward model
        Pfwd = getSimVER(Phi0, Peig, Fwd, sim, ap, ti, verbose)
#%% Step 2) Observe Forward Model (create vector of observations)
        bn = getObs(sim,cam,Lfwd,ti,Pfwd,makeplot,verbose)
        drnAll.append(bn)
#%% Step 3) fit constituent energies to our estimated vHat and reproject

        if overrides and overrides['fwdguess'][0] == 'maxwellian':
            Phi0z = maxwellian(Peig['Ek'],1e3,1e10)[0].ravel(order='F')
            Phi0r = outer(Phi0z, getpx(Fwd['x'],Wkm=1e3,X0=0,xs='gaussian'))
        elif overrides and overrides['fwdguess'][0] == 'true':
            Phi0r = Phi0.ravel(order='F')
            warn('** WARNING: feeding minimizer the true answer--testing only!! **')
        elif overrides and overrides['fwdguess'][0] == 'randn':
            randfact = absolute(normal(1,overrides['fwdguess'][1], size=Phi0.size))
            warn('** WARNING: feeding minizer true answer times {}'.format(randfact))
            Phi0r = randfact * Phi0.ravel(order='F')
        else: #normal case, no a priori
            Phi0r = zeros(Fwd['sx']*Peig['Mp'].shape[1]) #ones() is NOT appropriate -- must be tapered down for higher energy beams to keep physically plausible.

        Pfit,jfit,Tm,bfit = FitVER(Lfwd, bn, Phi0r, Peig, sim, cam,Fwd, ti, makeplot,verbose)
#%% collect results
        PhifitAll.append(jfit); bfitAll.append(bfit)
#%% plot results
        goPlot(sim,Fwd,cam,Lfwd,Tm,bn,bfit,Pfwd,Pfit,Peig,Phi0,
                     jfit,rawdata,ti,makeplot,progms,x1d,vlim,verbose)
        if animtime is not None:
            draw()
            pause(animtime)
        elif 'show' in makeplot:
            show()
        else:
            close('all')

#%% wrapup
    msg='{} done looping'.format(argv[0]); print(msg); #print(msg,file=stderr)

    png2multipage(progms,'.eps','.tif',descr=cmd,delete=False,verbose=verbose) #gif writing is not working yet

    analyseres(sim,cam,Fwd['x'],Fwd['xPixCorn'],
                   Phi0all,PhifitAll,drnAll,bfitAll,vlim,
                   x0true=None,E0true=None,
                   makeplot=makeplot,progms=progms,verbose=verbose)
#%% debug: save variables to MAT file
    if 'mat' in makeplot and progms is not None:
        from scipy.io import savemat
        cMatFN = join(progms,'comparePy.mat')
        try:
            print('saving to:',cMatFN)
            vd = {'drnP':bn,'LP':Lfwd,'vP':Pfwd,'vfitP':Pfit,#'vhatP':Phat['vART'],
                  'xPixCorn':Fwd['xPixCorn'],'zPixCorn':Fwd['zPixCorn']}
            savemat(cMatFN,oned_as='column',mdict=vd )
        except Exception as e:
            warn('failed to save to mat file.  {}'.format(e))

    msg ='{} program end'.format(argv[0]); print(msg); #print(msg,file=stderr)

    return Phi0all,PhifitAll #keep these, used for registration.py


def signal_handler(signal, frame):
    exit('\n *** Aborting program as per user pressed Ctrl+C ! \n')
#%% -----------------------------------------------------------
if __name__ == '__main__':
    from argparse import ArgumentParser

    signal(SIGINT, signal_handler) #allows friendly ctrl-c interrupt

    p = ArgumentParser(description='analyzes HST data and makes simulations')
    p.add_argument('infile',help='.xls filename with simulation parameters')
    p.add_argument('outdir',help='directory for output')
    p.add_argument('-v','--verbose',help='set debugging verbosity e.g. -v -vv -vvv',action='count',default=0)
    p.add_argument('--mat',help='save matlab .mat file of results',action="store_true")
    p.add_argument('--h5',help='save HDF5 .h5 file of fit results',action="store_true")
    p.add_argument('-m','--makeplot',help='list plots you want made',nargs='+',default=['']) #None gave type errors in generators & list comp.
    p.add_argument('-p','--showplot',help='show plots on screen',action="store_true")
    p.add_argument('-f','--frames',help='START STOP STEP of time indices',nargs=3,type=int)
    p.add_argument('--profile',help='Profile performance of program (development/debug only)',action='store_true')
    p.add_argument('-c','--cam',help='zero-indexed cameras to use (overrides XLS)',nargs='+',default=[None],type=int)
    p.add_argument('--cx',help='override cam positions (must specify all)',nargs='+',default=[None],type=float)
    p.add_argument('--influx',help='flux .h5 file to use (overrides XLS)')
    p.add_argument('--logplot',help='logarithmic axis scaling where appropriate',action='store_true')
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
    if ar.h5: makeplot.append('h5')
    if ar.mat: makeplot.append('mat')
    if ar.showplot: makeplot.append('show')
    if ar.logplot: makeplot.append('log')
    # these matplotlib imports MUST GO IN THIS ORDER
    import matplotlib as mpl
    #from mpl_toolkits.mplot3d import Axes3D #causes  TypeError: unhashable type: 'list'
    if in1d(makeplot,('png','eps')).any():
        print('using Agg backend: Visibly displayed plots are disabled!')
        mpl.use('Agg') #for fast PNG writing, but does NOT display at all!
    elif in1d(makeplot,'pdf').any():
        print('using PDF backend: Visibly displayed plots are disabled!')
        mpl.use('pdf') #for multipage PDF writing, but does NOT display at all!
    else:
        pass
        #mpl.use('Qt4Agg') # NOT FOR ANACONDA3
        #mpl.use('Tkagg') # possibly faster than qt4agg
    from matplotlib.pyplot import show,draw,pause
    print(('matplotlib backend / version: ' + mpl.get_backend() +'  ' + mpl.__version__  ))
#%%
    vlim = {'p':ar.vlim,'j':ar.jlim,'b':ar.blim}
#%%
    timeInds = ar.frames
    doProfile = ar.profile
#%% output directory
    progms = ar.outdir
    try:
        makedirs(progms)#, exist_ok=True) #python 2.7 doesn't have exist_ok
    except (OSError,TypeError) as e:
        pass


    with open(join(progms,'cmd.log'),'w') as f:
        f.write(' '.join(argv))
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
                 'overrides,progms,ar.x1d,vlim,ar.anim,' '.join(argv),ar.verbose)',proffn)
        pstats.Stats(proffn).sort_stats('time','cumulative').print_stats(50)
        #binpath = expanduser('~/profile/')
        #sysCall = [binpath + 'gprof2dot.py','-f','pstats',profFN,'|','dot','-Tpng','-o','output.png']
        #print(sysCall)
        #po = Popen(sysCall, stdout=PIPE, cwd=binpath, shell=False)
        #so,serr = po.communicate(timeout=1) #timeout is incompatible with Python 2.
        #print(so.decode('utf8'))

    else: #normal
        doSim(ParamFN,makeplot,timeInds,overrides,progms,ar.x1d, vlim,ar.anim,' '.join(argv), ar.verbose)

    if 'show' in makeplot:
        show()
