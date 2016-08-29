import logging
from numpy import empty,empty_like,isnan,sin,cos,radians,append,diff,ones,outer
import numpy as np #need this here
from scipy.sparse import csc_matrix
import h5py
from time import time
#
from .nans import nans
from .EllLineLength import EllLineLength

def getObs(sim,cam,L,tDataInd,ver):
    """
    real data: extract brightness vector from disk data
    simulation: create brightness from projection matrix and fwd model VER
    """
    if not sim.realdata and ver is None: #perhaps loading radar or allsky data
        return

    nCutPix = sim.ncutpix

    if sim.realdata:
        bn = empty(nCutPix * sim.nCamUsed,dtype=float,order='F') #FIXME assumes all cuts same length AND that cam 0 is used
        for C in cam:
            if not C.usecam: continue
            """
             remember that we put "d" in lexigraphical form,
             "d" is a column-major vector, such that if our 1D cut is N pixels,
             HST0 occupies d(0:N-1), HST1 occupies d(N:2N-1), ...
            """
            if C.keo.ndim==2: #more than 1 frame extracted
                thisCamPix = C.keo[:,tDataInd]
            elif C.keo.ndim==1:
                thisCamPix = C.keo
            else:
                raise ValueError('ndim==2 or 1 for stack of 1-D extracted cut')

            thisCamPix = mogrifyData(thisCamPix, C) #for clarity

            ''' here's where we assemble vector of observations'''
            bn[C.ind] = thisCamPix


    elif ver is not None: #or not np.any(np.isnan(v)): # no NaN in v # using synthetic data
        """ FIEFK """
 #       bp = L @ ver.ravel(order='F')
        bp = L.dot(ver.ravel(order='F'))
        assert bp.size == nCutPix * sim.nCamUsed

        bn = nans(bp.shape) #nans as a flag to check if something screwed up
        for C in cam:
            bn[C.ind] = mogrifyData(bp[C.ind],C)
#%% double check
   #assert np.any(np.isnan(drn)) == False # must be AFTER all drn are assigned, or you can get false positive errors!
    if isnan(bn).any():
        dumpFN = 'obsdump_t {}.h5'.format(tDataInd)
        logging.critical('NaN detected at tInd = {}   dumping variables to {}'.format(tDataInd,dumpFN))
        with h5py.File(str(dumpFN),'w',libver='latest') as fid:
            fid.create_dataset("/bn",data=bn)
            fid.create_dataset("/v",data=ver)
            fid.create_dataset("/L",data=L.todense(order='F'),compression='gzip')

    return bn

def makeCamFOVpixelEnds(Fwd,sim,cam,P):

    nCutPix = sim.ncutpix

    # when creating a new L, we always use all cameras, so ha['nCam'] is good here
    # order='F' not needed here because we don't reshape or ravel this variable
    xFOVpixelEnds = empty((nCutPix, sim.nCamUsed),dtype=float)
    zFOVpixelEnds = empty_like(xFOVpixelEnds)
#%% (1) define the two x,y points defining each pixel line of sight in the 2-D simulation grid
    for C in cam:
        '''LINE LENGTH
         here we have 2 vertices per angle instead of 3 (line instead of polygon)
         the minus sign on x makes the angle origin at local east
        '''
        if C.usecam:
            xFOVpixelEnds[:,C.name] = -(C.fovmaxlen * cos(radians(C.angle_deg))) + C.x_km
            zFOVpixelEnds[:,C.name] =  (C.fovmaxlen * sin(radians(C.angle_deg))) + C.alt_m/1000.

            C.xFOVpixelEnds = xFOVpixelEnds[:,C.name] #for plots.py
            C.zFOVpixelEnds = zFOVpixelEnds[:,C.name]

#%% (2) observational model auroral pixels
# now we make a matrix with the corner x,y coordinates of the auroral fwd
# model rectangular pixels

# x-coordinates of corners
#---------------
    #FIXME ASSUMES UNIFORM X-GRID
    Xpc = Fwd['x'] - 0.5 * sim.fwd_dxKM #shift left half a cell to get the corner
    Xpc = append(Xpc, Fwd['x'][-1] + 0.5 * sim.fwd_dxKM ) #last column by shifting half-cell to the right
    assert Xpc.size == Fwd['x'].size + 1
# y-coordinates of corners
#--------------
    if sim.useztranscar:
        Zpc = empty(Fwd['sz']+1, float)

        #FIXME would gradient() be better here?
        dz  = empty(Fwd['sz'],   float)
        dz[1:] = diff(Fwd['z'])
        dz[0] = dz[1] #FIXME assumes there was no grid jump (seems safe assumption)

        Zpc[:-1] = Fwd['z'] - 0.5*dz  # 0:Nz-1 values
        Zpc[-1] = Fwd['z'][-1] + 0.5*dz[-1] #FIXME assumes last two values have same spacing
    else:
        Zpc = Fwd['z'] - 0.5 * sim.fwd_dzKM #only for uniform z
        Zpc = append(Zpc, Fwd['z'][-1] + 0.5 * sim.fwd_dzKM) #last row

    assert Zpc.size == Fwd['z'].size + 1
    #we'll use these in doPlots and EllLineLength
    #no meshgrid anymore :-)
    Fwd['xPixCorn'] = Xpc
    Fwd['zPixCorn'] = Zpc

#%% (3) Compute intersection of Vol and FOV pixels (giving you "ell's")
# we say (for now) that ell=area of polygon intersection between FOV pixel and sky voxel
    tic = time()
    #used .values for future use of Numba
    L = EllLineLength(Fwd,xFOVpixelEnds,zFOVpixelEnds,
                      [c.x_km for c in cam if c.usecam],
                      [c.alt_m/1000. for c in cam if c.usecam],
                      nCutPix,sim,P['makeplot'],P['verbose'])
    print('computed L in {:0.1f}'.format(time()-tic) + ' seconds.')
    return L,Fwd,cam
#%%
def loadEll(sim,Fwd,cam,P):
    try:
      with h5py.File(str(sim.FwdLfn),'r',libver='latest') as fid:
        L = csc_matrix(fid['/L'])

        if Fwd is not None: #we're in main program
            if np.any(Fwd['x'] != fid['/Fwd/x']):#don't use .any() in case size is different
                raise ValueError('need to recompute L, as x-locations arent matched: loaded vs. commanded')
            if np.any(Fwd['z'] != fid['/Fwd/z']):  #don't use .any() in case size is different
                raise ValueError('need to recompute L, as z-locations arent matched: loaded vs. commanded')
        elif Fwd is None: #we're in another program
            Fwd = {'x':fid['/Fwd/x'].value,'z':fid['/Fwd/z'].value}
        else:
            raise TypeError("I dont understand what you're trying to do with Fwd, which should either be None or properly initialized")
        #{x,z}PixCorn must be assigned AFTER the if/elif/else
        Fwd['xPixCorn'] = fid['/Fwd/xPixCorn'].value
        Fwd['zPixCorn'] = fid['/Fwd/zPixCorn'].value

        if cam is not None:
            try:
                for i,C in enumerate(cam):
                    if C.usecam:
                        #cam[ci].angle_deg =  fid['/Obs/pixAngle'][:,i]
                        C.xFOVpixelEnds = fid['/Obs/xFOVpixelEnds'][:,i]
                        C.zFOVpixelEnds = fid['/Obs/zFOVpixelEnds'][:,i]
            except KeyError:
                logging.critical('could not load FOV ends, maybe this is an old Ell file')

        print('Loaded projection matrix L from {}'.format(sim.FwdLfn))

    except (FileNotFoundError,OSError) as e:
        logging.error('{} not found. Recomputing new Ell file. {}'.format(sim.FwdLfn,e))
        sim.loadfwdL = False
        L,Fwd,cam = getEll(sim,cam,Fwd,P)
    except AttributeError as e:
        logging.error('grid mismatch detected. use --ell command line option to save new Ell file. {}'.format(e))

    print('loadEll: Loaded "L,Fwd,Obs" data from: {}'.format(sim.FwdLfn))

    return L,Fwd,cam

def mogrifyData(data,cam):
    #steps should be in this order!
    data = data.astype(float)
    data *= cam.intens2dn #considers pixel area and camera amplifier gain
    data = cam.scaleintens(data) #camera cross-calibration

    data = cam.donoise(data)
    data = cam.dosmooth(data)
    data = cam.dolowerthres(data)

    data = cam.debias(data)
    data = cam.fixnegval(data)

    return data

def getEll(sim,cam,Fwd,P):

    if not sim.loadfwdL:
        if sim.nCamUsed != sim.useCamBool.size:
            logging.warning('To make a fresh L matrix, you must enable all HiST Cameras')

        L,Fwd,cam = makeCamFOVpixelEnds(Fwd,sim,cam,P)
    else:
        L,Fwd,cam = loadEll(sim,Fwd,cam,P)

    L = removeUnusedCamera(L,sim.useCamBool,sim.ncutpix)
    cam = definecamind(cam)

    assert 'x' in Fwd and Fwd['x'] is not None,'problem loading or computing grid'

    return L,Fwd,cam

def removeUnusedCamera(L,useCamBool,ncutpix):
    ''' remove unused cameras (rows of L) '''
    arow = ones(ncutpix,bool)
    grow = outer(arow,useCamBool).ravel(order='F')
    return L[grow,:]

def definecamind(cam):
    ''' store indices of b vector corresponding to each camera (in case some cameras not used) '''
    for i,C in enumerate(cam):
        if C.usecam:
            C.ind = slice(i*C.ncutpix, (i+1)*C.ncutpix)
    return cam
