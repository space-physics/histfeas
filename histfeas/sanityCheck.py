"""
sanity check for HiST parameters
Michael Hirsch
"""
from pathlib import Path
import logging
from pandas import read_excel
from shutil import copy2
from geopy.distance import vincenty
#
from .camclass import Cam
from .simclass import Sim


def getParams(XLSfn,overrides,makeplot,odir):
    XLSfn = Path(XLSfn).expanduser()
    if odir is not None:
        copy2(str(XLSfn),str(odir))
#%% read spreadsheet
    #paramSheets = ('Sim','Cameras','Arc')
    xl = read_excel(str(XLSfn),sheetname=None,index_col=0,header=0)
    sp = xl['Sim']
    cp = xl['Cameras']
#%% read arcs (if any)
    ap = {}; ntimeslice=None
    for s in xl:
        if s.startswith('Arc'):
            ap[s] = xl[s]
            if ntimeslice is not None and ap[s].shape[1]-1 != ntimeslice:
                raise ValueError('for now, all Arcs must have same number of times (columns)')
            ntimeslice=ap[s].shape[1]-1

    logging.info('# of observer time steps in spreadsheet: {}'.format(ntimeslice))
#%% class with parameters and function
    sim = Sim(sp,cp,ap,ntimeslice,overrides,makeplot,odir)
#%% grid setup
    Fwd = sim.setupFwdXZ(sp)
#%% setup cameras
    cam,cp = setupCam(sim,cp,Fwd['z'][-1])

    if sim.realdata:
        #find the x-coordinate (along B-perp) of each camera (can't do this inside camclass.py)
        cam[0].x_km = 0 # NOTE arbitrarily picks the first camera x=0km
        for C in cam[1:]:
            if C.usecam:
                C.x_km = vincenty((cam[0].lat,cam[0].lon),(C.lat,C.lon)).kilometers

    #store x,z in sim
    ellname=sim.getEllHash(sp,cp, [C.x_km for C in cam if C.usecam],[C.alt_m/1000. for C in cam if C.usecam])
    #will try to load this and compute if needed. Will be copied to output directory too.
    sim.FwdLfn = sim.rootdir/'precompute' / ellname

    # make the simulation time step match that of the fastest camera
    sim.kineticsec = min([C.kineticsec for C in cam if C.usecam])


    logging.info('fwd model voxels:\n'
          'B_perp: N={}   B_parallel: M={}'.format(Fwd['sx'],Fwd['sz']))
#%% init variables
    return ap,sim,cam,Fwd
###############################################

def setupCam(sim,cp,zmax):
    cam = []

    if sim.camxreq[0] is not None:
        logging.warning('overriding camera x-loc with {}'.format(sim.camxreq))
        for C,cx in zip(cp,sim.camxreq): #enumerate as in general camera 0 may not be used
            cp.at['Xkm',C] = cx
            cam.append(Cam(sim,cp[C], C, zmax))
    else:
        for C in cp:
            cam.append(Cam(sim,cp[C], C, zmax))

    assert len(cam)>0,'0 cams are configured, Nothing to do.'

    return cam,cp
