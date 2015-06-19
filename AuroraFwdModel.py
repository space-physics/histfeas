#!/usr/bin/env python3
"""
generates modeled aurora volume emission rate
Michael Hirsch
GPLv3+
"""
from __future__ import print_function, division
from numpy import isnan,zeros
from warnings import warn
#
from arcclass import getver
from transcararc import getColumnVER

def getSimVER(Phi0,Mp,Fwd,sim,ap,tInd,verbose):
    if not sim.realdata:
        Pfwd=zeros((Fwd['z'].size,Phi0.shape[1]),float) #TODO for non transcar cases
        for a in ap:
            #verify that all have same shape in B_\parallel
            assert (ap[a].loc['Zshape',0] == ap[a].loc['Zshape',:]).all()
            if ap[a].at['Zshape',0] == 'transcar':
                '''
                recall that phi(z,E) is 2-D matrix, and other matrices involved
                must be updated accordingly from original April 2014 try code
                '''
                #WARNING: The line below IS ****NOT**** += !
                return getColumnVER(sim.useztranscar, Mp['ztc'], Mp['Mp'], Phi0) 
            else:
                warn('these cases not yet updated for upsampled time!')
                Pfwd += getver(Fwd['x'], Fwd['z'], Mp, Phi0,
                          ap[a]['Wkm'], ap[a]['Hkm'],ap[a]['X0km'],
                          ap[a]['Z0km'],ap[a]['Xshape'], ap[a]['Zshape'],
                          ap[a]['Pnorm'])

        if isnan(Pfwd).any():
            return None
        else:
            return Pfwd
    else:
        return None #to ensure purity s.t. no unexpected ver is passed around when we're using real data
