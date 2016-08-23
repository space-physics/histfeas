#!/usr/bin/env python3
from numpy import outer, zeros_like,zeros,arange,repeat,nan, empty
#
from gridaurora.chapman import chapman_profile
from histutils.findnearest import find_nearest
from .transcararc import getpx


class Arc():

    def __init__(self,xl,texp):
        self.zshape = xl['Zshape']
        self.xshape = xl['Xshape']

        self.texp = texp

        for k in ('X0km','Pnorm','E0','Q0','Wbc','bl','bm','bh','Bm0','Bhf','Wkm'):
            try:
                self.__dict__[k] = getntimes(xl[k],texp.size)
            except KeyError:
                pass


def getntimes(sreq,N=None):
    """
    if 3 element req, see if it's a range spec.
    if 1 element, stretch to all time
    if empty elements, replace empties with nans
    """
    lreq = sreq.split(',') #preserves blank entries, numpy.fromstring doesn't preserve
    req = empty(len(lreq))
    for i,r in enumerate(lreq): # list of strings
        if len(r)==0: #empty string, passes 0
            req[i] = nan
        else:
            req[i] = float(r)
#%%
    if len(req)==3 and req[2] != 0 and abs(req[2]) < abs(req[1]):
        print('assuming .ini specifying value range')
        v = arange(req[0], req[1]+req[2],req[2])
    elif len(req) == 1 and N is not None: #replicate for all times
        v = repeat(req,N)
    else:
        v = req

    assert v.size >= 2,'must be at least 2 values to make 1 sim time step.'

    return v

def getver(x,z,Mp,Phi0, w,h,x0,z0,xshape,zshape,pmax):
    if zshape=='chapman':
        return ChapmanArc(w, h, x0,z0, x, z,xshape,pmax)
    elif zshape == 'rect':
        return RectArc(w,h, x0,z0, x,z, xshape,pmax)[0]
    elif zshape == 'transcar':
        return TypeError('use gettranscar arc function')
    elif zshape == 'zero': #zeros, already set to zero
        return zeros((z.size,x.size))
    else:
        raise ValueError('Unknown z-model: {}'.format(zshape))

def ChapmanArc(Wkm,H,X0,Z0,xKM,zKM,xshape, PC0=1):
    # chapman vert
    pz = PC0 * chapman_profile(Z0,zKM)
    # horizontal model
    px = getpx(xKM,Wkm,X0,xshape)
    # 2D model output
    ''' Make 2D Auroral Blob (original idea JLS)'''
    return outer(pz, px)

def RectArc(Wkm,Hkm,X0,Z0,xKM,zKM,xshape, PC0=1):

    #find lower and upper indices of rect. phantom
    PCind = (find_nearest( zKM, Z0-Hkm/2)[0], find_nearest(zKM, Z0+Hkm/2)[0] )
    #initialize vertical vector
    pz = zeros_like(zKM)
    pz[PCind[0]:PCind[1]+1] = PC0
    #horiz model
    px = getpx(xKM,Wkm,X0, xshape)

    pzx  = outer(pz, px)
    return pzx,pz,px
