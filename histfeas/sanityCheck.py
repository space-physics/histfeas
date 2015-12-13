"""
sanity check for HST simulation parameters
Michael Hirsch
GPL v3+

REQUIRES *** PANDAS 0.16 *** or newer for read_excel to work properly!
"""
from __future__ import division,absolute_import
from pathlib2 import Path
import logging
from pandas import read_excel
from shutil import copy2
from geopy.distance import vincenty
#
from .camclass import Cam
from .simclass import Sim


def getParams(XLSfn,overrides,makeplot,odir):
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
#%% ***** must be outside camclass ********
    nCutPix = cp.loc['nCutPix'].values
    if not (nCutPix == nCutPix[0]).all():
        raise ValueError('sanityCheck: all cameras must have same 1D cut length')
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
            C.x_km = vincenty((cam[0].lat,cam[0].lon),(C.lat,C.lon)).kilometers

    #store x,z in sim
    if odir and overrides and overrides['ell']:
        sim.FwdLfn = Path(odir) / sim.getEllHash(sp,cp,
                            [c.x_km for c in cam],[c.alt_m/1000. for c in cam])
    else:
        sim.FwdLfn = Path('precompute') / sim.getEllHash(sp,cp,
                            [c.x_km for c in cam],[c.alt_m/1000. for c in cam])

    # make the simulation time step match that of the fastest camera
    sim.kineticsec = min([C.kineticsec for C,u in zip(cam,sim.useCamBool) if u])


    logging.info('fwd model voxels:\n'
          'B_perp: N={}   B_parallel: M={}'.format(Fwd['sx'],Fwd['sz']))
#%% init variables
    return ap,sim,cam,Fwd
###############################################

def setupCam(sim,cp,zmax):
    cam = []

    if sim.camxreq[0] is not None:
        logging.warning('overriding camera x-loc with {}'.format(sim.camxreq))
        for i,(c,cx) in enumerate(zip(cp,sim.camxreq)):
            if sim.useCamBool[i]:
                cp.iat['Xkm',c] = cx
                cam.append(Cam(sim,cp[c], c, zmax))
    else:
        for i,c in enumerate(cp):
            if sim.useCamBool[i]:
                cam.append(Cam(sim,cp[c], c, zmax))

    if len(cam)==0:
        raise ValueError('0 cams are configured, Nothing to do.')
    return cam,cp
