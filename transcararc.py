from numpy import (asfortranarray,atleast_3d, exp,sinc,pi,zeros, outer,
                   isnan,log,logspace,where,empty,allclose)
import h5py
from scipy.interpolate import interp1d
from warnings import warn
from datetime import timedelta
#
from transcarutils.readTranscar import getTranscar
from transcarutils.eFluxGen import fluxgen
from histutils.findnearest import find_nearest


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
        warn('* cubic interpolating Transcar altitude, use caution that VER peaks arent missed...')
        fver = interp1d(zTranscar, Peig, axis=0, kind='cubic')
        Tm = asfortranarray(fver(zKM))

    return Tm.dot(Phi0)

def getMp(sim,zKM,makeplot,dbglvl):
#%% read from transcar sim
    Peigen,EKpcolor = getTranscar(sim, dbglvl)[:2]
    try:
        Ek = Peigen.columns.values
        zTranscar = Peigen.index.values
    except AttributeError as e:
        warn('*** getMp: it appears there was trouble earlier on with getTranscar, aborting. {}'.format(e))
        return None
#%% clip to Hist requested altitudes
    if not allclose(zKM,zTranscar):
        warn('** getMp: warning, attempting to trim altitude grid, this may not be successful due to floating point error')
        goodAltInd = (zKM[0] < zTranscar) &  (zTranscar < zKM[-1])
        Peig = asfortranarray(Peigen.values[goodAltInd,:])
    else:
        Peig = asfortranarray(Peigen.values)
#%% repack, with optional downsample
    if sim.downsampleEnergy:
        print('** downsampling in energy **')
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

def getPhi0(sim,ap,xKM,Ek,makeplots,dbglvl):
    #%% get flux
    if not sim.realdata:
        if sim.Jfwdh5 is not None:
            print('Loading sim. input diff. number flux from', sim.Jfwdh5)
            with h5py.File(sim.Jfwdh5,'r',libver='latest') as f:
                Phi0 = asfortranarray(atleast_3d(f['/phiInit']))
        else:
            Phi0 = empty((Ek.size,xKM.size,sim.nArc),order='F')
            pz = fluxgen(Ek, ap, dbglvl)[0]
    #%% horizontal modulation
            px = getpx(xKM,ap.loc['Wkm'].values.astype(float),
                           ap.loc['X0km'].values.astype(float),
                           ap.loc['Xshape'].values)
            for i in range(sim.nArc):
                Phi0[...,i]  = outer(pz[:,i], px[i,:])

        assert xKM.size == Phi0.shape[1]
    #%% obtain simulation time steps from spreadsheet
        tsim = []
        for i,ts in enumerate(ap.loc['tReqOffsetSec'].values.astype(float)):
            tsim.append(sim.transcarutc + timedelta(seconds=ts))
    else:
        Phi0 = None
    return Phi0,tsim

def getpx(xKM,Wkm,X0,xs='gaussian'):
    px = zeros((X0.size,xKM.size),order='F') #since numpy 2-D array naturally iterates over rows
#%%
    px[xs=='gaussian',:] = exp(-((xKM-X0[xs=='gaussian',None])/Wkm[xs=='gaussian',None])**2) #(original idea JLS)
#%%
    ir = where(xs=='rect')[0]
    for i in ir:
        #find leftmost and rightmost indices of rect. phantom
        PGind = (find_nearest(xKM, X0[i]-Wkm[i]/2)[0], find_nearest(xKM, X0[i]+Wkm[i]/2)[0] )
        px[i, PGind[0]:PGind[1]+1] = 1.
#%%
    px[xs=='sinc2',:] = sinc( pi*(xKM - X0[xs=='sinc2',None])/Wkm[xs=='sinc2',None] )**2
#%%
    px[xs=='none',find_nearest(xKM,X0[xs=='none'])] = 1.
#    jn = where(xs=='none')
#    for j in jn:
#        px[j,find_nearest(xKM,X0[j])[0]] = 1.

    return px
