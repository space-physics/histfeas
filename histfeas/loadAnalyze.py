#!/usr/bin/env python3
"""
To generate inputs for this program, run main_hist.py with
-m h5
option
"""

from __future__ import division,absolute_import
import h5py
from os.path import join,expanduser,splitext,isfile
from numpy import asarray,diff
from warnings import warn
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='poster', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'})
#
from analysehst import analyseres
from sanityCheck import getParams
from plotsnew import plotB, plotJ, plotVER
from observeVolume import definecamind

vlim={'b':(None,None),'j':(None,None),'p':[None]*6}

def runtest(h5list,xlsfn,overrides,makeplot,verbose=0):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]; Pest=[]; Pfwd=[]

    nt = len(h5list)
    if nt==0:
        warn('no HDF5 files found, ending.')
        return

    for h5 in h5list:
        with h5py.File(h5,'r',libver='latest') as f:
            x  = f['/pfwd/x'].value #same for all in directory
            xp = f['/pfwd/xp'].value

            Phifwd.append(f['/phifwd/phi'].value)
            Phidict.append({'x':f['/phiest/phi'].value,
                            'EK':f['/phiest/Ek'].value})

            Pest.append(f['/pest/p'].value)
            Pfwd.append(f['/pfwd/p'].value)
            z = f['/pest/z'].value
            zp = f['/pest/zp'].value

            dhat.append(f['/best/bfit'].value)
            drn.append(f['/best/braw'].value)
#%%
    stem,ext = splitext(h5)  #all in same directory

    Phifwd = asarray(Phifwd).transpose(1,2,0)

    if xlsfn is None:
        warn('No XLSX parameter file found')
        return

    ap,sim,cam,Fwd = getParams(xlsfn,overrides,makeplot,progms=None,verbose=verbose)
    cam = definecamind(cam,sim.nCutPix)

    for a in ap:
        #TODO assumes all are same distance apart
        x0true = ap[a].loc['X0km',:][:-1] + 0.5*diff(ap[a].loc['X0km',:])
        E0true = ap[a].loc['E0',:][:-1]   + 0.5*diff(ap[a].loc['E0',:])

        analyseres(None,None,
                   x, xp, Phifwd, Phidict, drn, dhat,
                   vlim, x0true,E0true,makeplot,
                   progms=ext, verbose=verbose)

#%%
    for ti in range(nt):
        if 'fwd' in makeplot:
            plotB(drn[ti],sim.realdata,cam,vlim['b'],ti,makeplot,'$br',None,verbose)

            plotJ(sim,Phifwd[...,ti],x,xp,Phidict[ti]['EK'],None,
                  vlim['j'],vlim['p'][:2],ti,makeplot,'phifwd',
                  '$\phi_{top}$ fwd diff. number flux',
                  None,None,verbose)

            plotVER(sim,Pfwd[ti],x,xp,z,zp,vlim['p'],ti,makeplot,'pfwd',
                    '$P$ fwd volume emission rate',
                    None,None,verbose)

        if 'optim' in makeplot:
            plotJ(sim,Phidict[ti]['x'],x,xp,Phidict[ti]['EK'],None,
                  vlim['j'],vlim['p'][:2],ti,makeplot,'phiest',
                  '$\hat{\phi}_{top}$ estimated diff. number flux',
                  None,None,verbose)

            plotVER(sim,Pest[ti],x,xp,z,zp,vlim['p'],ti,makeplot,'pest',
                    '$\hat{P}$ estimated volume emission rate',
                    None,None,verbose)

if __name__ == '__main__':
    from glob import glob
    from os.path import dirname
    #
    from argparse import ArgumentParser
    p = ArgumentParser(description="load HiST output and plot/analyse")
    p.add_argument('h5path',help='path containing dump.h5 outputs (Hist output)')
    p.add_argument('-m','--makeplot',help='plots to make',default=[])
    p = p.parse_args()

    h5path = expanduser(p.h5path)

    if isfile(h5path):
        h5list = [h5path]
        xlsfn = glob(join(dirname(h5path),'*.xlsx'))
    else:
        h5list = glob(join(h5path,'dump_*.h5'))
        h5list.sort()
        xlsfn = glob(join(h5path,'*.xlsx'))

    if len(xlsfn) == 0: xlsfn=None

    runtest(h5list,xlsfn[0],None,p.makeplot)

    show()