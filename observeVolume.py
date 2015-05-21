from numpy import empty,empty_like,s_,isnan,sin,cos,radians,append,diff,ones,outer
import numpy as np
from scipy.sparse import csc_matrix
from EllLineLength import EllLineLength
import h5py
from warnings import warn
from time import time
from nans import nans
#from pdb import set_trace

def getObs(sim,cam,L,tDataInd,ver,makePlots,dbglvl):
    """
    real data: extract brightness vector from disk data
    simulation: create brightness from projection matrix and fwd model VER
    """

    nCutPix = sim.nCutPix

    if sim.realdata:
        bn = empty(nCutPix * sim.nCamUsed,dtype=float,order='F') #FIXME assumes all cuts same length AND that cam 0 is used
        for c in cam:
            cInd = cam[c].ind
            """
             remember that we put "d" in lexigraphical form,
             "d" is a column-major vector, such that if our 1D cut is N pixels,
             HST0 occupies d(0:N-1), HST1 occupies d(N:2N-1), ...
            """
            thisCamPix  = cam[c].keo[:,tDataInd]

            thisCamPix = mogrifyData(thisCamPix, cam[c]) #for clarity

            ''' here's where we assemble vector of observations'''
            bn[cInd] = thisCamPix


    elif ver is not None: #or not np.any(np.isnan(v)): # no NaN in v # using synthetic data
        """ FIEFK """
        bp = L.dot(ver.ravel(order='F'))
        assert bp.size == nCutPix * sim.nCamUsed

        bn = nans(bp.shape) #nans as a flag to check if something screwed up
        for c in cam:
            cInd = cam[c].ind
            bn[cInd] = mogrifyData(bp[cInd],cam[c])

    else: #skip VER processing
      print('skipping VER generation and pixel projection due to None in VER')
      bn = None

#%% double check
   #assert np.any(np.isnan(drn)) == False # must be AFTER all drn are assigned, or you can get false positive errors!
    if bn is not None and isnan(bn).any():
        dumpFN = 'obsdump_t {}.h5'.format(tDataInd)
        warn('NaN detected at tInd = {}   dumping variables to {}'.format(tDataInd,dumpFN))
        with h5py.File(dumpFN,'w',libver='latest') as fid:
            fid.create_dataset("/bn",data=bn)
            fid.create_dataset("/v",data=ver)
            fid.create_dataset("/L",data=L.todense(order='F'),compression='gzip')

    return bn

def makeCamFOVpixelEnds(Fwd,sim,cam,makePlots,dbglvl):

    nCutPix = sim.nCutPix

    # when creating a new L, we always use all cameras, so ha['nCam'] is good here
    # order='F' not needed here because we don't reshape or ravel this variable
    xFOVpixelEnds = empty((nCutPix, sim.nCamUsed),dtype=float)
    zFOVpixelEnds = empty_like(xFOVpixelEnds)
#%% (1) define the three x,y points defining each 2D pixel cone
    for c in cam:#.keys():
        '''LINE LENGTH
         here we have 2 vertices per angle instead of 3 (line instead of polygon)
         the minus sign on x makes the angle origin at local east
        '''
        xFOVpixelEnds[:,int(c)] = -(cam[c].fovmaxlen * cos(radians(cam[c].angle_deg) ) +
                   cam[c].x_km )
        zFOVpixelEnds[:,int(c)] =  (cam[c].fovmaxlen * sin(radians(cam[c].angle_deg) ) +
                   cam[c].z_km)

        cam[c].xFOVpixelEnds = xFOVpixelEnds[:,int(c)] #for plots.py
        cam[c].zFOVpixelEnds = zFOVpixelEnds[:,int(c)]

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
        Zpc = empty(Fwd['sz']+1, dtype=float)

        #FIXME would gradient() be better here?
        dz  = empty(Fwd['sz'],   dtype=float)
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
                            sim.allCamXkm,
                            sim.allCamZkm,
                            sim.savefwdL,nCutPix,sim,makePlots,dbglvl)
    print('computed L in {:0.1f}'.format(time()-tic) + ' seconds.')
    return L,Fwd,cam
#%%
def loadEll(Fwd,cam,EllFN,verbose):
    try:
      with h5py.File(EllFN,'r',libver='latest') as fid:
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
                for i,c in enumerate(cam):
                    #cam[ci].angle_deg =  fid['/Obs/pixAngle'][:,i]
                    cam[c].xFOVpixelEnds = fid['/Obs/xFOVpixelEnds'][:,i]
                    cam[c].zFOVpixelEnds = fid['/Obs/zFOVpixelEnds'][:,i]
            except KeyError:
                warn('could not load FOV ends, maybe this is an old Ell file')

    except (IOError) as e: #python 2.7 doesn't have FileNotFoundError
      raise IOError('*** loadEll: {} not found.\nuse --ell command line option to save new Ell file. {}'.format(EllFN,e))
    except AttributeError as e:
        raise AttributeError('grid mismatch detected. use --ell command line option to save new Ell file. {}'.format(e))

    if verbose>0: print('loadEll: Loaded "L,Fwd,Obs" data from:', EllFN)

    return L,Fwd,cam

def mogrifyData(data,cam):
    #steps should be in this order!
    data *= cam.intens2dn #considers pixel area and camera amplifier gain
    data = cam.scaleintens(data) #camera cross-calibration

    data = cam.donoise(data)
    #print('noise std. deviation cam ',cam.name,'=',cam.std)
    data = cam.dosmooth(data)
    data = cam.dolowerthres(data)

    data = cam.debias(data)
    data = cam.fixnegval(data)

    return data

def getEll(sim,cam,Fwd,makePlots,dbglvl):

    if not sim.loadfwdL:
        if sim.nCamUsed != sim.useCamBool.size:
            raise ValueError('To make a fresh L matrix, you must enable ALL cameras all(useThisCam == 1) ')

        L,Fwd,cam = makeCamFOVpixelEnds(Fwd,sim,cam,makePlots,dbglvl)
    else:
        L,Fwd,cam = loadEll(Fwd,cam,sim.FwdLfn,dbglvl)

    L,cam = removeUnusedCamera(L,sim.useCamBool,sim.nCutPix,cam)

    return L,Fwd,cam

def removeUnusedCamera(L,useCamBool,nCutPix,cam):
    ''' remove unused cameras (rows of L) '''
    arow = ones(nCutPix).astype(bool)
    grow = outer(arow,useCamBool).ravel(order='F')
    L = L[grow,:]
    ''' store indices of b vector corresponding to each camera (in case some cameras not used) '''
    for i,c in enumerate(cam):
        cam[c].ind = s_[ i*nCutPix : (i+1)*nCutPix ]

    return L,cam
