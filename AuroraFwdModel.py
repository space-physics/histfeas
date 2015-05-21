#!/usr/bin/env python3
"""
generates modeled aurora volume emission rate
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division
from numpy import isnan
#
from arcclass import getver

def getSimVER(Phi0,Mp,Fwd,sim,ap,tInd,verbose):
    if not sim.realdata:
        Pfwd = getver(Fwd['x'], Fwd['z'], Mp, Phi0,
                      ap['Wkm'], ap['Hkm'],ap['X0km'],
                      ap['Z0km'],ap['Xshape'], ap['Zshape'],
                      sim.useztranscar,ap['Pnorm'])

        if isnan(Pfwd).any():
            Pfwd = None
    else:
        Pfwd = None #to ensure purity s.t. no unexpected ver is passed around when we're using real data

    return Pfwd