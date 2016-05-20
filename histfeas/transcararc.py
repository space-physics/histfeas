from numpy import (asfortranarray,atleast_3d, exp,sinc,pi,zeros, outer,in1d,
                   isnan,log,logspace,arange,allclose,diff,atleast_1d,isfinite,empty,nan)
import h5py
from scipy.interpolate import interp1d
import logging
from xarray import DataArray
#
from gridaurora.eFluxGen import fluxgen
from gridaurora.arcexcite import getTranscar
from histutils.findnearest import find_nearest

def getColumnVER(zgrid,zTranscar,Peig,Phi0):
    assert Phi0.shape[0] == Peig.shape[1]
    assert zTranscar.shape[0] == Peig.shape[0]

    if zgrid: #using original transcar z-locations
        Tm = Peig
    else:
        raise NotImplementedError('this interpolation was rarely used so disabled it.')
#        warn('* cubic interpolating Transcar altitude, use caution that VER peaks arent missed...')
#        fver = interp1d(zTranscar, Peig, axis=0, kind='cubic')
#        Tm = asfortranarray(fver(zKM))

    return Tm @ Phi0

def getMp(sim,cam,zKM,makeplot):
    if not in1d(('fwd','optim'),makeplot).any():
        return
#%% read from transcar sim
    Peigen,EKpcolor = getTranscar(sim,cam[0].alt_m/1000.,90-cam[0].Bincl)[:2]
    assert isinstance(Peigen,DataArray),'Did not get DataArray from getTranscar, aborting.'
    Ek = Peigen.energy_ev.values
    zTranscar = Peigen.alt_km.values
#%% clip to Hist requested altitudes
    if not allclose(zKM,zTranscar):
        logging.warning('attempting to trim altitude grid, this may not be successful due to floating point error')
        goodAltInd = (zKM[0] < zTranscar) &  (zTranscar < zKM[-1])
        Peig = asfortranarray(Peigen.values[goodAltInd,:])
    else:
        Peig = asfortranarray(Peigen.values)
#%% repack, with optional downsample
    if sim.downsampleEnergy:
        logging.warning('** downsampling in energy **')
        Ek,EKpcolor,Peigen = downsampleEnergy(Ek,EKpcolor,Peig, sim.downsampleEnergy)
    #FIXME: just use a DataFrame!
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
        logging.warning('should these NaNs be set to zero?')
    return Ek2,EKpcolor2,Mp2

def getPhi0(sim,arc,xKM,Ek,makeplots):
    #%% get flux
    if not sim.realdata:
        if sim.Jfwdh5 is not None:
            print('Loading sim. input diff. number flux from {}'.format(sim.Jfwdh5))
            with h5py.File(str(sim.Jfwdh5),'r',libver='latest') as f:
                Phi0 = asfortranarray(atleast_3d(f['/phiInit']))
        else:
            Phi0 = assemblePhi0(sim,arc,Ek,xKM)
        assert xKM.size == Phi0.shape[1]
    else:
        Phi0 = None
    return Phi0

def assemblePhi0(sim,arc,Ek,xKM):
    Phi0 = zeros((Ek.size,xKM.size,sim.nTimeSlice),order='F') #NOT empty, since we sum to build it!

    for a in arc: #iterate over arcs, using superposition
#%% upsample to sim time steps
        E0,Q0,Wbc,bl,bm,bh,Bm,Bhf, Wkm,X0 = upsampletime(arc[a],sim)

        if arc[a].zshape == 'transcar':
            phiz = fluxgen(Ek, E0,Q0,Wbc,bl,bm,bh,Bm,Bhf)[0] # Nenergy x Ntime
        elif arc[a].zshape == 'flat':
            phiz = zeros((Ek.size, sim.nTimeSlice)) #zeros not empty or nan
            for i,e in enumerate(E0):
                try:
                    phiz[Ek<=e,i] = Q0[i]  # Nenergy x Ntime_sim
                except ValueError:
                    pass
        elif arc[a].zshape == 'impulse':
            phiz = zeros((Ek.size, sim.nTimeSlice)) #zeros not empty or nan
            for i,e in enumerate(E0):
                try:
                    phiz[find_nearest(Ek,e)[0],i] = Q0[i]  # Nenergy x Ntime_sim
                except ValueError:
                    pass
        else:
            raise NotImplementedError('unknown zshape = {}'.format(arc.zshape))
#%% horizontal modulation
        phix = getpx(xKM,Wkm,X0,arc[a].xshape)

        if arc[a].zshape == 'transcar':
            for i in range(sim.nTimeSlice):
                #unsmeared in time
                phi0sim = zeros((Ek.size,xKM.size),order='F') #NOT empty, since we're summing!
                for j in range(sim.timestepsperexp):
                    phi0sim += outer(phiz[:,i*sim.timestepsperexp+j],
                                     phix[i*sim.timestepsperexp+j,:])
                Phi0[...,i] += phi0sim
        elif arc[a].zshape in ('impulse','flat'):
            phix[~isfinite(phix)] = 0.
            for i in range(sim.nTimeSlice):
                Phi0[...,i] += outer(phiz[:,i],phix[i,:])
        else:
            raise NotImplementedError

    return Phi0

def upsampletime(arc,sim):
    #%% obtain observation time steps from spreadsheet (for now, equal to kinetic time)
    if abs(sim.kineticsec - diff(arc.texp).mean()) > 1e-3:
        logging.error('exposure time not matching spreadsheet arc time step')
    # make simulation time, also defined as seconds since Transcar tReq
    dtsim =sim.kineticsec/sim.timestepsperexp
    tsim = arange(arc.texp[0],arc.texp[-1],dtsim)

# FUTURE
#    #tsim is a finer time step than texp, the camera exposure
#    tsim = empty(texp.size*sim.timestepsperexp,dtype=datetime)
#    tsimstep = timedelta(seconds=sim.kineticsec/sim.timestepsperexp)
#    for i,t in enumerate(texp):
#        #sim time steps (for future, in case spreadsheet steps != to exposure time (kinetic time))
#        for j in range(i*sim.timestepsperexp, (i+1)*sim.timestepsperexp):
#            tsim[j] = sim.transcarutc + j*tsimstep

    #probably could be done with lambda

    f = interp1d(arc.texp, arc.E0);     E0  = f(tsim)
    if (~isfinite(E0)).any():
        E0 = arc.E0[:-1]

    f = interp1d(arc.texp, arc.Q0);     Q0  = f(tsim)
    if (~isfinite(Q0)).any():
        Q0 = arc.Q0[:-1]

    f = interp1d(arc.texp, arc.Wbc);  Wbc = f(tsim)
    f = interp1d(arc.texp, arc.bl);   bl  = f(tsim)
    f = interp1d(arc.texp, arc.bm);   bm  = f(tsim)
    f = interp1d(arc.texp, arc.bh);   bh  = f(tsim)
    f = interp1d(arc.texp, arc.Bm0);  Bm  = f(tsim)
    f = interp1d(arc.texp, arc.Bhf);  Bhf = f(tsim)

    f = interp1d(arc.texp, arc.Wkm);    Wkm = f(tsim)
    if (~isfinite(Wkm)).any():
        Wkm = arc.Wkm[:-1]

    f = interp1d(arc.texp, arc.X0km);   X0  = f(tsim)
    if (~isfinite(X0)).any():
        X0 = arc.X0km[:-1]

    logging.info('new E0 upsamp [eV]: {}'.format(E0))

    return E0,Q0,Wbc,bl,bm,bh,Bm,Bhf, Wkm, X0

def getpx(xKM,Wkm,X0,xs):
    assert isinstance(xs,str)
    X0=atleast_1d(X0); Wkm=atleast_1d(Wkm)
    px = zeros((X0.size,xKM.size),order='F') #since numpy 2-D array naturally iterates over rows
#%%
    if xs =='gaussian':
        px = exp(-((xKM-X0[:,None])/Wkm[:,None])**2) #(original idea JLS)
#%%
    elif xs =='rect':
    #ir = where(xs=='rect')[0]
    #for i in ir:
        for i in range(X0.size):
            #find leftmost and rightmost indices of rect. phantom
            PGind = (find_nearest(xKM, X0[i]-Wkm[i]/2)[0],
                     find_nearest(xKM, X0[i]+Wkm[i]/2)[0] )
            px[i, PGind[0]:PGind[1]+1] = 1.
#%%
    elif xs == 'sinc2':
        px= sinc( pi*(xKM - X0[:,None])/Wkm[:,None] )**2
#%%
    else:
        px[:,find_nearest(xKM,X0)] = 1.
#    jn = where(xs=='none')
#    for j in jn:
#        px[j,find_nearest(xKM,X0[j])[0]] = 1.

    return px
