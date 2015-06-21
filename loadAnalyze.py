#!/usr/bin/env python3
from __future__ import division,absolute_import
import h5py
from os.path import join,expanduser,splitext
from numpy import asarray
from warnings import warn
from matplotlib.pyplot import show
#
from analysehst import analyseres
from sanityCheck import getParams
from plotsnew import plotB, plotJ
from observeVolume import definecamind

vlim={'b':(None,None),'j':(None,None),'p':(None,None)}

def runtest(h5list,xlsfn,overrides,makeplot,verbose=0):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]
    
    nt = len(h5list)
    
    for h5 in h5list:
        with h5py.File(h5,'r',libver='latest') as f:
            x  = f['/pfwd/x'].value #same for all in directory
            xp = f['/pfwd/xp'].value
            
            Phifwd.append(f['/phifwd/phi'].value)
            Phidict.append({'x':f['/phiest/phi'].value, 
                            'EK':f['/phiest/Ek'].value})

            dhat.append(f['/best/bfit'].value)
            drn.append(f['/best/braw'].value)
#%%  
    stem,ext = splitext(h5)  #all in same directory
    
    Phifwd = asarray(Phifwd).transpose(1,2,0)
        
    analyseres(None,None,
               x, xp, Phifwd, Phidict, drn, dhat,
               vlim=vlim, makeplot=[None],
               progms=ext, verbose=verbose)

#%%
    if xlsfn is None: 
        warn('No XLSX parameter file found')
        return
        
    ap,sim,cam,Fwd = getParams(xlsfn,overrides,makeplot,progms=None,verbose=verbose)
    cam = definecamind(cam,sim.nCutPix)
 
    for ti in range(nt):
        plotB(drn[ti],sim.realdata,cam,vlim['b'],ti,makeplot,'$br',None,verbose)
        
        plotJ(sim,Phidict[ti]['x'],x,xp,Phidict[ti]['EK'],None,
          vlim['j'],vlim['p'][:2],ti,makeplot,'phiest',
          '$\hat{\phi}_{top}$ estimated diff. number flux',
          None,None,verbose)
    
if __name__ == '__main__':
    from glob import glob    
    #
    from argparse import ArgumentParser
    p = ArgumentParser(description="load HiST output and plot/analyse")
    p.add_argument('h5path',help='path containing dump.h5 outputs (Hist output)')
    p = p.parse_args()
    
    h5list = glob(expanduser(join(p.h5path,'dump_*.h5')))
    h5list.sort()
    
    xlsfn = glob(expanduser(join(p.h5path,'*.xlsx')))
    if len(xlsfn) == 0: xlsfn=[None]
    
    runtest(h5list,xlsfn[0],overrides=None,makeplot=('fwd','optim'))
    
    show()