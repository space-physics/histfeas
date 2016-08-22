#!/usr/bin/env python3
from numpy import outer, zeros_like,zeros,fromstring,arange
#
from gridaurora.chapman import chapman_profile
from histutils.findnearest import find_nearest
from .transcararc import getpx


class Arc():

    def __init__(self,xl,texp):
        self.zshape = xl['Zshape']
        self.xshape = xl['Xshape']

        self.texp = texp

        self.X0km = getntimes(xl['x0km'])

        self.Pnorm = fromstring(xl['Pnorm'],sep=',')

        self.E0 = fromstring(xl['E0'],sep=',')
        self.Q0 = fromstring(xl['Q0'],sep=',')
        self.Wbc = fromstring(xl['Wbc'],sep=',')
        self.bl = fromstring(xl['bl'],sep=',')
        self.bm = fromstring(xl['bm'],sep=',')
        self.bh = fromstring(xl['bh'],sep=',')
        self.Bm0 = fromstring(xl['Bm0'],sep=',')
        self.Bhf = fromstring(xl['Bhf'],sep=',')


        self.Wkm = fromstring(xl['Wkm'],sep=',')

def getntimes(req):
    """
    if 3 element req, see if it's a range spec. Otherwise, take literally
    """
    req = fromstring(req,sep=',')
    if len(req)==3 and abs(req[2]) < abs(req[1]):
        print('assuming .ini specifying value range')
        v = arange(req[0], req[1]+req[2],req[2])
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
        raise TypeError('Unknown model type: {}'.format(zshape))

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
