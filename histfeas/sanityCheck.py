"""
sanity check for HST simulation parameters
Michael Hirsch
GPL v3+

REQUIRES *** PANDAS 0.16 *** or newer for read_excel to work properly!
"""
from pandas import read_excel
from warnings import warn
from shutil import copy2
#
from .camclass import Cam
from .simclass import Sim

def getParams(XLSfn,overrides,makeplot,progms,verbose):
    if progms is not None:
        copy2(XLSfn,progms)
#%% read spreadsheet
    #paramSheets = ('Sim','Cameras','Arc')
    xl = read_excel(XLSfn,sheetname=None,index_col=0,header=0)
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

    print('# of observer time steps in spreadsheet: {}'.format(ntimeslice))
#%% ***** must be outside camclass ********
    nCutPix = cp.loc['nCutPix'].values
    if not (nCutPix == nCutPix[0]).all():
        raise ValueError('sanityCheck: all cameras must have same 1D cut length')
#%% class with parameters and function
    sim = Sim(sp,cp,ap,ntimeslice,overrides,makeplot,progms,verbose)
#%% grid setup
    Fwd = sim.setupFwdXZ(sp)
#%% setup cameras
    cam,cp = setupCam(sim,cp,Fwd['z'][-1],verbose)

    print('fwd model voxels:\n'
          'B_perp: N={}   B_parallel: M={}'.format(Fwd['sx'],Fwd['sz']))
#%% init variables
    return ap,sim,cam,Fwd
###############################################

def setupCam(sim,cp,zmax,dbglvl):
    cam = {}

    if sim.camxreq[0] is not None:
        warn('overriding camera x-loc with {}'.format(sim.camxreq))
        for i,(c,cx) in enumerate(zip(cp,sim.camxreq)):
            if sim.useCamBool[i]:
                cp.iat['Xkm',c] = cx
                cam[c] = Cam(sim,cp[c], c, zmax,dbglvl)
    else:
        for i,c in enumerate(cp):
            if sim.useCamBool[i]:
                cam[c] = Cam(sim,cp[c], c, zmax,dbglvl)

    if len(cam)==0:
        raise ValueError('0 cams are configured, Nothing to do.')
    return cam,cp
