"""
sanity check for HiST parameters
Michael Hirsch
"""
from . import Path
import logging
from six.moves.configparser import ConfigParser
from shutil import copy2
from geopy.distance import vincenty
#
from .camclass import Cam
from .simclass import Sim
from .arcclass import Arc

def getParams(inifn,overrides,makeplot,odir):
    inifn = Path(inifn).expanduser()
    if odir is not None:
        copy2(str(inifn),str(odir))
#%% read spreadsheet
    #paramSheets = ('Sim','Cameras','Arc')
    try:
        xl = ConfigParser(allow_no_value=True, inline_comment_prefixes=('#'), strict=True)
    except TypeError: #py27
        xl = ConfigParser(allow_no_value=True)
    xl.read(str(inifn))
#%% read arcs (if any)
    arc,ntimeslice = setupArc(xl)
    logging.info('# of observer time steps in spreadsheet: {}'.format(ntimeslice))
#%% class with parameters and function
    sim = Sim(xl,arc,ntimeslice,overrides,makeplot,odir)
#%% grid setup
    Fwd = sim.setupFwdXZ(xl)
#%% setup cameras
    cam = setupCam(sim,xl['cam'],Fwd['z'][-1],makeplot)

    if sim.realdata:
        #find the x-coordinate (along B-perp) of each camera (can't do this inside camclass.py)
        cam[0].x_km = 0 # NOTE arbitrarily picks the first camera x=0km
        for C in cam[1:]:
            if C.usecam:
                C.x_km = vincenty((cam[0].lat,cam[0].lon),(C.lat,C.lon)).kilometers

    #store x,z in sim
    ellname=sim.getEllHash(xl, [C.x_km for C in cam],[C.alt_m/1000. for C in cam])
    #will try to load this and compute if needed. Will be copied to output directory too.
    sim.FwdLfn = sim.rootdir/'precompute' / ellname

    # make the simulation time step match that of the fastest camera
    sim.kineticsec = min([C.kineticsec for C in cam if C.usecam])


    logging.info('fwd model voxels:\n'
          'B_perp: N={}   B_parallel: M={}'.format(Fwd['sx'],Fwd['sz']))
#%% init variables
    return arc,sim,cam,Fwd
###############################################

def setupArc(xl):
    arc = {}
    ntimeslice=None

    print(xl)

    for s in xl.sections(): # for py27
        if s.startswith('arc'):
            if ntimeslice is not None and getntimes(xl[s]) != ntimeslice:
                raise ValueError('for now, all Arcs must have same number of times (columns)')
            ntimeslice = getntimes(xl[s]) # last time is blended with 2nd to last time

            arc[s] = Arc(xl[s])

    return arc, ntimeslice

def getntimes(arc):
    return len(arc['tsec'].split(','))-1

def setupCam(sim,cp,zmax,makeplot):
    cam = []

    if sim.camxreq[0] is not None:
        logging.warning('overriding camera x-loc with {}'.format(sim.camxreq))
        for C,cx in zip(sim.camnames,sim.camxreq): #enumerate as in general camera 0 may not be used
            cam.append(Cam(sim, cp, C, zmax,makeplot))
    else:
        for C in sim.camnames:
            cam.append(Cam(sim, cp, C, zmax,makeplot))

    assert len(cam)>0,'0 cams are configured, Nothing to do.'

    return cam
