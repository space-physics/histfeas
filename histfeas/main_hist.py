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
from sys import argv
from os.path import join
from os import makedirs
from numpy import absolute,zeros,outer
from numpy.random import normal
from warnings import warn
from matplotlib.pyplot import close,draw,pause,show
#
from gridaurora.eFluxGen import maxwellian
from histutils.imageconv import png2multipage
from histutils.simulFrame import getSimulData
from histutils.get1Dcut import get1Dcut #we need cam.angle_deg for plotting
from .sanityCheck import getParams
from .AuroraFwdModel import getSimVER
from .transcararc import getMp,getPhi0,getpx #calls matplotlib
from .observeVolume import getEll,getObs #calls matplotlib
from .FitVER import FitVERopt as FitVER #calls matplotlib
from .plotsnew import goPlot #calls matplotlib
from .analysehst import analyseres
#import logging
#logging.basicConfig(filename='hist.log',filemode='w',level=logging.DEBUG)


def doSim(ParamFN,makeplot,timeInds,overrides,progms,x1d,vlim,animtime, cmd,verbose):

    #%% output directory
    try:
        makedirs(progms)#, exist_ok=True) #python 2.7 doesn't have exist_ok
    except (OSError,TypeError) as e:
        pass

    with open(join(progms,'cmd.log'),'w') as f:
        f.write(' '.join(argv))
#%% Step 0) load data
    ap,sim,cam,Fwd = getParams(ParamFN, overrides,makeplot,progms,verbose)
#%% setup loop
    if sim.realdata:
        cam,rawdata,sim = getSimulData(sim,cam,makeplot,progms,verbose)
    else: #simulation
        rawdata = None
        if sim.raymap == 'astrometry':
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
    else:
        Phi0all = None
    if verbose>0: print('timeInds: {}'.format(timeInds))
#%%start looping for each time slice in keogram (just once if simulated)
    for ti in timeInds:
        print('entering time {}'.format(ti))
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
        try:
            if overrides['fwdguess'][0] == 'maxwellian':
                Phi0z = maxwellian(Peig['Ek'],1e3,1e10)[0].ravel(order='F')
                Phi0r = outer(Phi0z, getpx(Fwd['x'],Wkm=1e3,X0=0,xs='gaussian'))
            elif overrides['fwdguess'][0] == 'true':
                Phi0r = Phi0.ravel(order='F')
                warn('** WARNING: feeding minimizer the true answer--testing only!! **')
            elif overrides['fwdguess'][0] == 'randn':
                randfact = absolute(normal(1,overrides['fwdguess'][1], size=Phi0.size))
                warn('** WARNING: feeding minizer true answer times {}'.format(randfact))
                Phi0r = randfact * Phi0.ravel(order='F')
            else: #normal case, no a priori
                Phi0r = zeros(Fwd['sx']*Peig['Mp'].shape[1]) #ones() is NOT appropriate -- must be tapered down for higher energy beams to keep physically plausible.
        except (KeyError,TypeError):
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
