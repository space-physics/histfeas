#!/usr/bin/env python3
"""
generates modeled aurora volume emission rate
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division
from numpy import isnan
#
from arcclass import Arc

def getSimVER(Phi0,Mp,Fwd,sim,ap,tInd,dbglvl):
    if not sim.realdata:
        arc = Arc(sim,ap)

        Pfwd = arc.getver(Fwd['x'], Fwd['z'], Mp, Phi0)

        if isnan(Pfwd).any():
            Pfwd = None
    else:
        Pfwd = None #to ensure purity s.t. no unexpected ver is passed around when we're using real data
        arc = None

    return Pfwd,arc