"""
sanity check for HST simulation parameters
Michael Hirsch
GPL v3+
"""
from pandas import read_excel
from camclass import Cam
from simclass import Sim
from xlrd.biffh import XLRDError

def getParams(XLSfn,overrides,savedump,progms,dbglvl):
#%% read spreadsheet
    paramSheets = ('Sim','Cameras','Arcs')
    sp = read_excel(XLSfn,paramSheets[0],index_col=0,header=0)
    cp = read_excel(XLSfn,paramSheets[1],index_col=0,header=0)
    try:
        ap = read_excel(XLSfn,paramSheets[2],index_col=0,header=0)
    except XLRDError:
        ap = None
#%% ***** must be outside camclass ********
    nCutPix = cp.loc['nCutPix'].values
    if not (nCutPix == nCutPix[0]).all():
        exit('*** all cameras must have same 1D cut length')
#%% class with parameters and function
    sim = Sim(sp,cp,ap,overrides,progms,dbglvl)
#%% grid setup
    Fwd = sim.setupFwdXZ(sp)
#%% setup cameras
    cam,cp = setupCam(Cam,sim,cp,Fwd['z'][-1],dbglvl)

    print('fwd model voxels:\n')
    print('B_perp: N=',Fwd['sx'],'   B_parallel: M=',Fwd['sz'])
#%% init variables
    return ap,sim,cam,Fwd
###############################################

def setupCam(Cam,sim,cp,zmax,dbglvl):
    cam = {}

    if sim.camxreq[0] is not None:
        print('* overriding camera x-loc with ',sim.camxreq)
        for i,ci,cx in zip(sim.useCamInd,sim.useCamInd.astype(str),sim.camxreq):
            cp.ix['Xkm',i] = cx
            cam[ci] = Cam(sim,cp.ix[:,i], ci, zmax,dbglvl)
    else:
        for i,ci in zip(sim.useCamInd,sim.useCamInd.astype(str)):
            cam[ci] = Cam(sim,cp.ix[:,i], ci, zmax,dbglvl)
    return cam,cp
