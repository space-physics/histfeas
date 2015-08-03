"""
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division, absolute_import
from numpy import exp,outer, zeros_like,zeros
from pybashutils.findnearest import find_nearest

from transcararc import getpx

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
    pz = PC0 * exp(.5*(1-(zKM-Z0)/H - exp((Z0-zKM)/H)))
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