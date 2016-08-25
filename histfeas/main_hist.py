#!/usr/bin/env python
"""
 Michael Hirsch

 Procedure:
 0) ingest simulation parameters from XLS file, sanity check, load data
 1) Make an auroral phantom (if not real data, a Simulation, and optionally for a first guess for reconstruction)
 2) create projection matrix
 3) fit differential number flux to reconstructed VER
 4) Test reconstruction & fit by reprojecting to synthetic image brightness

example: (fwd model only)
python main_hist.py in/2cam_trans.xlsx /dev/shm/rev_trans2/ -m fwd png --vlim -0.5 3.5 90 350 1e9 1e10 --jlim 1e5 5e5 --blim 0 1e4 -f 0 120 20
"""
import logging
from sys import argv
from numpy import absolute,zeros,outer
from numpy.random import normal
from time import time
from matplotlib.pyplot import close,draw,pause,show
#
from gridaurora.eFluxGen import maxwellian
from histutils.simulFrame import getSimulData
from histutils.get1Dcut import get1Dcut #we need cam.angle_deg for plotting
#
from . import getParams
from .AuroraFwdModel import getSimVER
from .transcararc import getMp,getPhi0,getpx #calls matplotlib
from .observeVolume import getEll,getObs #calls matplotlib
from .FitVER import FitVERopt as FitVER #calls matplotlib
from .plotsnew import goPlot #calls matplotlib

def doSim(P):
    print('')
#%% housekeeping
    tic = time()
    logging.basicConfig(level=30-P['verbose']*10)
    #%% output directory
    P['outdir'].mkdir(parents=True,exist_ok=True)

    (P['outdir']/'cmd.log').write_text(' '.join(argv)) #store command for future log

#%% Step 0) load data
    arc,sim,cam,Fwd,P = getParams(P)
#%% setup loop
    if sim.realdata:
        # can load enormous amount of data into rawdata, Ncam x Nframe x Ny x Nx (verify dim order?)
        cam,rawdata,sim = getSimulData(sim,cam,P['outdir'],P['verbose'])
    else: #simulation
        rawdata = None
        if sim.raymap == 'astrometry':
            cam = get1Dcut(cam,P['outdir'],P['verbose'])
    timeInds = sim.maketind(P['timeinds'])
#%% Step 1) get projection matrix
    Lfwd,Fwd,cam = getEll(sim,cam,Fwd,P)
#%% load eigenprofiles from Transcar
    Peig = getMp(sim,cam,Fwd['z'],P['makeplot'])
#%% synthetic diff. num flux
    Phi0all = getPhi0(sim,arc,Fwd['x'],Peig['Ek'], P['makeplot']) # Nenergy x Nx x Ntime
    print('{:.1f} sec to prepare for HiSTfeas loop'.format(time()-tic))
#%%start looping for each time slice in keogram (just once if simulated)
    for ti in timeInds:
        logging.info('entering time {}'.format(ti))
        if sim.realdata:
            Phi0 = None; Pfwd = None
        else: #sim
            """
            we need to integrate in time over the relevant time slices
            the .sum(axis=2) does the integration/smearing in time
            """
            Phi0 = Phi0all[...,ti] # Nenergy x Nx
#%% Step 1) Forward model
        Pfwd = getSimVER(Phi0, Peig, Fwd, sim, arc, ti) # Nz x Nx
#%% Step 2) Observe Forward Model (create vector of observations)
        bn = getObs(sim,cam,Lfwd,ti,Pfwd) # Ncam*Npixel (1D vector)
#%% Step 3) fit constituent energies to our estimated vHat and reproject
        Phi0r = initPhi(Phi0,Peig,Fwd,P['overrides'])
        Pfit,jfit,Tm,bfit = FitVER(Lfwd, bn, Phi0r, Peig, sim, cam,Fwd, ti, P)
#%% plot results
        goPlot(sim,Fwd,cam,Lfwd,Tm,bn,bfit,Pfwd,Pfit,Peig,Phi0,
                     jfit,rawdata,ti,P)
        if 'animtime' in P and P['animtime']:
            draw()
            pause(P['animtime'])
        elif 'show' in P['makeplot']:
            show()
        else:
            close('all')

#%% wrapup
    msg ='{} program end'.format(argv[0]); print(msg); #print(msg,file=stderr)

def initPhi(Phi0,Peig,Fwd,overrides):
    try:
        if overrides['fwdguess'][0] == 'maxwellian':
            Phi0z = maxwellian(Peig['Ek'],1e3,1e10)[0].ravel(order='F')
            return outer(Phi0z, getpx(Fwd['x'],Wkm=1e3,X0=0,xs='gaussian'))
        elif overrides['fwdguess'][0] == 'true':
            logging.error('Feeding minimizer the true answer--testing only')
            return Phi0.ravel(order='F')
        elif overrides['fwdguess'][0] == 'randn':
            randfact = absolute(normal(1,overrides['fwdguess'][1], size=Phi0.size))
            logging.error('feeding minizer true answer times {}'.format(randfact))
            return randfact * Phi0.ravel(order='F')
        else: #normal case, no a priori
            return zeros(Fwd['sx']*Peig['Mp'].shape[1]) #ones() is NOT appropriate -- must be tapered down for higher energy beams to keep physically plausible.
    except (KeyError,TypeError):
        try:
            return zeros(Fwd['sx']*Peig['Mp'].shape[1]) #ones() is NOT appropriate -- must be tapered down for higher energy beams to keep physically plausible.
        except TypeError: # no fwd or optim
            pass
