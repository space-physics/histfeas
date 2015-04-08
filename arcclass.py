"""
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division
from numpy import exp,outer, zeros_like,zeros,isnan
from findnearest import find_nearest

from transcararc import getColumnVER,getpx

class Arc:
    def __init__(self,sim,ap):
        self.xshape = ap['Xshape'].lower()
        self.zshape = ap['Zshape'].lower()
        self.h = ap['Hkm']

        if isnan(ap['Wkm']):
            exit('*** you must specify W0 for all used arcs')

        if isnan(ap['X0km']):
            exit('*** you must specify X0 for all used arcs')

        if isnan(ap['Pnorm']):
            exit('*** you must specify P0 for all used arcs')


        self.w = ap['Wkm']
        self.x0 = ap['X0km']
        self.pmax = ap['Pnorm']
        self.z0 = ap['Z0km']
        self.xLim = sim.fwd_xlim

        self.zgrid = sim.useztranscar


    def getver(self, x, z, Mp, Phi0):
        if self.zshape=='chapman':
            return ChapmanArc(self.w, self.h, self.x0,self.z0,
                               x, z,
                               self.xshape,
                               self.pmax)[0]
        elif self.zshape == 'rect':
            return RectArc(self.w, self.h, self.x0,self.z0,
                                  x,z,
                                  self.xshape,
                                  self.pmax)[0]
        elif self.zshape == 'transcar':
            '''
            recall that phi(z,E) is 2-D matrix, and other matrices involved
            must be updated accordingly from original April 2014 try code
            '''
            return getColumnVER(self.zgrid, Mp['ztc'], Mp['Mp'], Phi0, z)

        elif self.zshape == 'zero': #zeros, already set to zero
            return zeros((z.size,x.size))
        else:
            exit('*** Unknown model type ' + self.zshape)

def ChapmanArc(Wkm,H,X0,Z0,xKM,zKM,xshape, PC0=1):
    # chapman vert
        pz = PC0 * exp(.5*(1-(zKM-Z0)/H - exp((Z0-zKM)/H)))
    # horizontal model
        px = getpx(xKM,Wkm,X0,xshape)
    # 2D model output
        ''' Make 2D Auroral Blob (original idea JLS)'''
        pzx  = outer(pz, px)
        return pzx,pz,px

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