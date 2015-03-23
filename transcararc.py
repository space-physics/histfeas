from numpy import (asfortranarray,atleast_3d, exp,sinc,pi,zeros_like, outer,
                   isnan,log,logspace,nan)
from sys import path
import h5py
from scipy.interpolate import interp1d
#
path.append('../transcar-utils')
from readTranscar import getTranscar
from eFluxGen import fluxgen
from findnearest import find_nearest

def getTranscarMp(sim,makeplot,dbglvl):
 #%% get VER/flux
    Peigen, EKpcolor = getTranscar(sim, dbglvl)[:2]

    #return Mp,zTranscar,Ek,EKpcolor
    return asfortranarray(Peigen.values), Peigen.index.values, Peigen.columns.values, EKpcolor

def getColumnVER(zgrid,zTranscar,Peig,Phi0,zKM):
    assert Phi0.shape[0] == Peig.shape[1]
    assert zTranscar.shape[0] == Peig.shape[0]

    if zgrid: #using original transcar z-locations
        Tm = Peig
    else:
        print('* cubic interpolating Transcar altitude, use caution that VER peaks arent missed...')
        fver = interp1d(zTranscar, Peig, axis=0, kind='cubic')
        Tm = asfortranarray(fver(zKM))

    return Tm.dot(Phi0)

def getMp(sim,zKM,makeplot,dbglvl):
#%% read from transcar sim
    Peigen,EKpcolor = getTranscar(sim, dbglvl)[:2]
    try:
        Ek = Peigen.columns.values; zTranscar = Peigen.index.values
    except AttributeError as e:
        exit('getMp: it appears there was trouble earlier on with getTranscar, aborting. ' + str(e))
#%% clip to Hist resqueted alittudes
    goodAltInd = (zTranscar <= zKM[-1]) #TODO clips top only
    Peig = asfortranarray(Peigen.values[goodAltInd,:])
#%% repack, with optional downsample
    if sim.downsampleEnergy:
        Ek,EKpcolor,Peigen = downsampleEnergy(Ek,EKpcolor,Peig, sim.downsampleEnergy)
    return {'Mp':Peig,'ztc':zTranscar,'Ek':Ek,'EKpcolor':EKpcolor}

def downsampleEnergy(Ek,EKpcolor,Mp,downsamp):
    """ we know original points are logspaced.
    1) make new Ek2 axis, with 1/downsamp as many log-spaced points
    2) fill new Phi0_downsamp values with sums of adjacent Ek bins
    """
    nEK = Ek.size//downsamp #floor

    Ek2 = logspace(log(Ek[0]), log(Ek[-1]),
                num=nEK, endpoint=True, base=exp(1))
    EKpcolor2 = logspace(log(EKpcolor[0]), log(EKpcolor[-1]),
                num=nEK+1, endpoint=True, base=exp(1))

    fp = interp1d(log(Ek), log(Mp), kind='linear',axis=1)
    Mp2 = exp(fp(log(Ek2)))
    if isnan(Mp2).any():
        print('** downsampleEnergy: should these NaNs be set to zero?')
    return Ek2,EKpcolor2,Mp2

def getPhi0(sim,ap,PhiFN,xKM,Ek,minev,makeplots,dbglvl):
    #%% get flux
    if not sim.realdata:
        if PhiFN is not None:
            if dbglvl>0:
                print('Loading sim. input diff. number flux from', PhiFN)
            with h5py.File(PhiFN,'r',libver='latest') as f:
                Phi0 = asfortranarray(atleast_3d(f['/phiInit']))
            Phi0[Ek<minev,:,:] = 0 #this is the FORTRAN 3-D axis order
        else:
            pz = fluxgen(Ek,ap['E0'],ap['Q0'],ap['Wbc'],ap['Bm'],ap['Bhf'],ap['bl'],
                               ap['bm'],ap['bh'],'fluxgen' in makeplots,dbglvl)[0]
    #%% horizontal modulation
            px = getpx(xKM,ap['Wkm'],ap['X0km'],
                       ap['Xshape'],ap['Pnorm'])
            Phi0  = outer(pz, px)
            Phi0[Ek<minev,:] = 0
    #%% manual zeroing of low fluxes
        assert xKM.size == Phi0.shape[1]
    else:
        Phi0 = None
    return Phi0

def getpx(xKM,Wkm,X0,BperpShape='gaussian',P0=1):
    if isnan(X0) or isnan(Wkm) or isnan(P0):
        print('*** getpx: you must define X0, W0 and P0')

    if BperpShape == 'gaussian':
        px = exp(-((xKM-X0)/Wkm)**2) #(original idea JLS)
    elif BperpShape == 'rect':
        px = zeros_like(xKM);
        #find leftmost and rightmost indices of rect. phantom
        PGind = (find_nearest(xKM, X0-Wkm/2)[0], find_nearest(xKM, X0+Wkm/2)[0] )

        px[PGind[0]:PGind[1]+1] = P0
    elif BperpShape == 'sinc2':
        px = sinc( pi*(xKM - X0)/Wkm )**2
    elif BperpShape == 'none': #no smear
        px = zeros_like(xKM)
        px[find_nearest(xKM,X0)[0]] = 1
    else:
        print('*** Unknown B_perp shape ' + str(BperpShape))
        return nan

    return px
