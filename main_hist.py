#!/usr/bin/env python3
"""
 Michael Hirsch
 May 2014
 GPL v3+

 Procedure:
 0) ingest simulation parameters from XLS file, sanity check, load data
 1) Make an auroral phantom (if not real data, a Simulation, and optionally for a first guess for reconstruction)
 2) create projection matrix
 3) reconstruction of VER
 4) fit differential number flux to reconstructed VER
 5) Test reconstruction & fit by reprojecting to synthetic image brightness

"""
from __future__ import division,print_function
from signal import signal,SIGINT #for Ctrl C
from os.path import expanduser, join
from os import makedirs
from numpy import absolute,zeros,asarray,in1d,arange
from numpy.random import normal
import h5py
from datetime import timedelta
#
from sanityCheck import getParams

def doSim(ParamFN,savedump,makeplot,datadump,timeInds,overrides,progms,x1d,vlim,animtime, dbglvl):
    # local -- these were put here so that matplotlib backend autoselect could happen first
    from AuroraFwdModel import getSimVER
    from transcararc import getMp,getPhi0 #calls matplotlib
    from observeVolume import getEll,getObs #calls matplotlib
    from FitVER import FitVERopt as FitVER #calls matplotlib
#    from simulFrame import getSimulData
    from plotsnew import goPlot
    from analysehst import analyseres
#%% Step 0) load data
    ap,sim,cam,Fwd = getParams(ParamFN, overrides,savedump,makeplot,progms,dbglvl)
#%% setup loop
    if sim.realdata:
        cam,rawdata,sim = getSimulData(sim,cam,makeplot,progms,dbglvl)
        nTimeSlice = cam['0'].keo.shape[1] #FIXME assumes equal num. of time slices list(cam)[0]
    else: #simulation
        rawdata = None
        if sim.raymap == 'astrometry': #and any('b' in m[:2] for m in makeplot):
            from get1Dcut import get1Dcut #we need cam.angle_deg for plotting
            cam = get1Dcut(cam,makeplot,progms,dbglvl)
        nTimeSlice = ap.shape[1]
    timeInds = sim.maketind(timeInds,nTimeSlice)
#%% Step 1) get projection matrix
    Lfwd,Fwd,cam = getEll(sim,cam,Fwd,makeplot,dbglvl)
#%% preallocation
    jfitAll = []; drnAll = []; bfitAll=[];  jAll= []; PfitAll = []
#%% load eigenprofiles from Transcar
    Peig = getMp(sim,Fwd['z'],makeplot,dbglvl)
#%% synthetic diff. num flux
    if not sim.realdata:
        try:
            Phi0all,tsim = getPhi0(sim,ap,Fwd['x'],Peig['Ek'], makeplot,dbglvl)
        except TypeError as e:
            print('*** trouble with forward model. Exiting sim.   {}'.format(e))
            return
#%%start looping for each time slice in keogram (just once if simulated)
    for ti in timeInds:
        if sim.realdata:
            Phi0 = None; Pfwd = None; arc = None
        else: #sim
            """
            we need to integrate in time over the relevant time slices
            the .sum(axis=2) does the integration/smearing in time
            """
            #FIXME qualify our selected times from the command line, so we don't overlap in time
            cti = [tsim.index(t) for t in tsim if t>=tsim[ti]  and t<(tsim[ti]+timedelta(seconds=sim.kineticSec) ) ]
            Phi0 = Phi0all[...,cti].sum(axis=2)
#%% Step 1) Forward model
        Pfwd,arc = getSimVER(Phi0, Peig, Fwd, sim, ap.iloc[:,ti], ti, dbglvl)
#%% Step 2) Observe Forward Model (create vector of observations)
        bn = getObs(sim,cam,Lfwd,ti,Pfwd,makeplot,dbglvl)
        drnAll.append(bn)
#%% Step 3) Reconstruction (not used in long time)
        #vhat,vr0 = ReconVER(L,bn,sim,ap,cp,Fwd,ti,makeplot,dbglvl)
#%% Step 4) fit constituent energies to our estimated vHat and reproject
        if overrides['fwdguess'][0] == 'true':
            Phi0r = Phi0.ravel(order='F')
            print('****************************************************************')
            print('** WARNING: feeding minimizer the true answer--testing only!! **')
            print('****************************************************************')
        elif overrides['fwdguess'][0] == 'randn':
            randfact = absolute(normal(1,overrides['fwdguess'][1], size=Phi0.size))
            print('** WARNING: feeding minizer true answer times {}'.format(randfact))
            Phi0r = randfact * Phi0.ravel(order='F')
        else: #normal case, no a priori
            Phi0r = zeros(Fwd['sx']*Peig['Mp'].shape[1]) #ones() is NOT appropriate -- must be tapered down for higher energy beams to keep physically plausible.

        Pfit,jfit,Tm,bfit = FitVER(Lfwd, bn, Phi0r, Peig, sim, cam,Fwd, makeplot,dbglvl)
#%% collect results
        jfitAll.append(jfit); bfitAll.append(bfit)
#%% plot results
        goPlot(ParamFN,sim,arc,Fwd,cam,Lfwd,Tm,bn,bfit,Pfwd,Pfit,Phi0,
                     jfit,rawdata,ti,makeplot,progms,x1d,vlim,dbglvl)
        if animtime is not None:
            plt.draw()
            plt.pause(animtime)
        elif 'show' in makeplot:
            plt.show()
        else:
            plt.close('all')

        if 'h5' in savedump:
            PfitAll.append(Pfit)
#%%
    print('done looping')
    analyseres(sim,Fwd['x'],Fwd['xPixCorn'],cam,
                   Phi0all,jfitAll,drnAll,bfitAll,vlim,makeplot,progms)
#%% debug: save variables to MAT file
    if 'mat' in savedump:
        from scipy.io import savemat
        cMatFN = ''.join((progms,'/comparePy','.mat'))
        try:
            print('saving to:',cMatFN)
            vd = {'drnP':bn,'LP':Lfwd,'vP':Pfwd,'vfitP':Pfit,#'vhatP':Phat['vART'],
                  'xPixCorn':Fwd['xPixCorn'],'zPixCorn':Fwd['zPixCorn']}
            savemat(cMatFN,oned_as='column',mdict=vd )
        except Exception as e:
            print('failed to save to mat file.  {}'.format(e))
    if 'h5' in savedump:
        ch5fn = ''.join((progms,'/fitted','.h5'))
        jAll = asarray(jAll).transpose(axes=(1,2,0)) #asarray puts pages first
        print('saving to:', ch5fn)
        with h5py.File(ch5fn,'w',libver='latest') as f:
            f["/Jflux"]=jAll
            f["/Vfit"]=PfitAll

#%% =======================================================================
def signal_handler(signal, frame):
    print('\n *** Aborting program as per user pressed Ctrl+C ! \n')
    exit(0)
#%% -----------------------------------------------------------
if __name__ == '__main__':
    from argparse import ArgumentParser
    from sys import argv

    signal(SIGINT, signal_handler) #allows friendly ctrl-c interrupt

    p = ArgumentParser(description='analyzes HST data and makes simulations')
    p.add_argument('infile',help='.xls filename with simulation parameters',type=str)
    p.add_argument('outdir',help='directory for output',type=str,nargs='?',default='out')
    p.add_argument('-d','--debug',help='set debugging verbosity',default=0,type=float)
    p.add_argument('--mat',help='save matlab .mat file of results',action="store_true")
    p.add_argument('--h5',help='save HDF5 .h5 file of fit results',action="store_true")
    p.add_argument('--dump',help='dump debugging data to .h5 at key points of program',action="store_true")
    p.add_argument('-m','--makeplot',help='list plots you want made',nargs='+',default=[''],type=str) #None gave type errors in generators & list comp.
    p.add_argument('-p','--showplot',help='show plots on screen',action="store_true")
    p.add_argument('-f','--frames',help='START STOP STEP of time indices',nargs=3,default=None,type=int)
    p.add_argument('--profile',help='Profile performance of program (development/debug only)',action='store_true')
    p.add_argument('-c','--cam',help='zero-indexed cameras to use (overrides XLS)',nargs='+',default=[None],type=int)
    p.add_argument('--cx',help='override cam positions (must specify all)',nargs='+',default=[None],type=float)
    p.add_argument('--influx',help='flux .h5 file to use (overrides XLS)',default=None,type=str)
    p.add_argument('--logplot',help='logarithmic axis scaling where appropriate',action='store_true')
    p.add_argument('--x1d',help='required location [km] of x for 1-D flux plots',default=None,type=float)
#    p.add_argument('--saveall',help='saves "all" variables to disk, useful for registration case tracking',action='store_true')
    p.add_argument('-a','--anim',help='animate plots (crudely)',type=float,default=None)
    p.add_argument('--filter',help='optical filter choices: bg3   none',type=str,default=None)
    p.add_argument('--minev',help='minimum beam energy to include in fwd model',type=float,default=None)
    p.add_argument('--fwdguess',help='feed minimizer fwd answer. true | randn stddev |',type=str,nargs='+',default=[None])
    p.add_argument('--fitm',help='override fit (minimization) method',type=str,default=None)
    p.add_argument('--vlim',help='xmin xmax zmin zmax pmin pmax   limits for VER plots',type=float,nargs=6,default=(None,None,None,None,None,None))
    p.add_argument('--jlim',help='MIN MAX flux limits for diff num flux plots',type=float,nargs=2,default=(None,None))
    p.add_argument('--blim',help='MIN MAX flux limits for brightness plots',type=float,nargs=2,default=(None,None))
    p.add_argument('-L','--ell',help='force recomputation of sparse matrix L',action='store_true')
    ar = p.parse_args()

    ParamFN = expanduser(ar.infile)
    savedump=[None]
    if ar.mat: savedump.append('mat')
    if ar.h5: savedump.append('h5')
#%% plot setup
    makeplot = ar.makeplot
    if ar.showplot: makeplot.append('show')
    if ar.logplot: makeplot.append('log')
    # these matplotlib imports MUST GO IN THIS ORDER
    import matplotlib as mpl
    #from mpl_toolkits.mplot3d import Axes3D #causes  TypeError: unhashable type: 'list'
    if in1d(makeplot,('png','eps')).any():
        print('** using Agg backend: Visibly displayed plots are disabled!')
        mpl.use('Agg') #for fast PNG writing, but does NOT display at all!
    elif in1d(makeplot,'pdf').any():
        print('** using PDF backend: Visibly displayed plots are disabled!')
        mpl.use('pdf') #for multipage PDF writing, but does NOT display at all!
    else:
        pass
        #mpl.use('Qt4Agg') # NOT FOR ANACONDA3
        #mpl.use('Tkagg') # possibly faster than qt4agg
    import matplotlib.pyplot as plt
    print(('matplotlib backend / version: ' + mpl.get_backend() +'  ' + mpl.__version__  ))
#%%
    vlim = {'p':ar.vlim,'j':ar.jlim,'b':ar.blim}
#%%
    timeInds = ar.frames
    doProfile = ar.profile
    datadump = ar.dump
#%% output directory
    progms = ar.outdir
    try:
        makedirs(progms)#, exist_ok=True) #python 2.7 doesn't have exist_ok
    except (OSError,TypeError) as e:
        pass #for python 2.7


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
        import cProfile
        from profileRunHSTsim0 import goCprofile
        profFN = 'hstprofstats.pstats'
        print('saving profile results to ' + profFN)
        cProfile.run('doSim(ParamFN,savedump,makeplot,datadump,timeInds,overrides,progms,ar.x1d,vlim,ar.anim,dbglvl=ar.debug)',profFN)

        #binpath = expanduser('~/profile/')
        #sysCall = [binpath + 'gprof2dot.py','-f','pstats',profFN,'|','dot','-Tpng','-o','output.png']
        #print(sysCall)
        #po = Popen(sysCall, stdout=PIPE, cwd=binpath, shell=False)
        #so,serr = po.communicate(timeout=1) #timeout is incompatible with Python 2.
        #print(so.decode('utf8'))
        goCprofile(profFN)
    else: #normal
        doSim(ParamFN,savedump,makeplot,datadump,timeInds,
                                      overrides,progms,ar.x1d, vlim,ar.anim,dbglvl=ar.debug)

#    if ar.saveall:
#        with h5py.File(progms + 'validate.h5',libver='latest') as fid:
#            fid.create_dataset('/fitp/residual',data=fitpAll['/optimresidual'])
#            fid.create_dataset("/L",data=L,        compression="gzip")
    #exit(progms)

    #if saveplots and makePlots:
   #     print('converting to TIFF')
   #     subprocess.call(['convert','out/Jadj_t*.png','-compress','zip','out/JadjSim.tiff'], shell=False)
