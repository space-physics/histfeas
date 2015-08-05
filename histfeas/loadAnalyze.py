#!/usr/bin/env python3
"""
To generate inputs for this program, run main_hist.py with
-m h5
option
"""

from __future__ import division,absolute_import
import h5py
from os import makedirs
from os.path import join,expanduser,splitext,isfile,dirname
from tempfile import gettempdir
from numpy import asarray,diff
from warnings import warn
from matplotlib.pyplot import show
from glob import glob
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='poster', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'})
#
try:
    from .analysehst import analyseres
    from .sanityCheck import getParams
    from .plotsnew import plotB, plotJ, plotVER,plotBcompare
    from .observeVolume import definecamind
except:
    from analysehst import analyseres
    from sanityCheck import getParams
    from plotsnew import plotB, plotJ, plotVER,plotBcompare
    from observeVolume import definecamind

vlim={'b':(None,None),'j':(None,None),'p':[None]*6}

def readresults(h5list,xlsfn,overrides,makeplot,verbose=0):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]; Pest=[]; Pfwd=[]
    tInd = [];

    nt = len(h5list)
    if nt==0:
        warn('no HDF5 files found, ending.')
        return

    for h5 in h5list:
        stem = splitext(h5)[0]

        tInd.append(int(stem[-3:])) #NOTE assumes last 3 digits are time ind
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
    stem,ext = splitext(h5)  #all in same directory, left here for clarity
    progms = join(gettempdir(),dirname(h5))
    try:
        makedirs(progms)
    except:
        pass

    Phifwd = asarray(Phifwd).transpose(1,2,0)

    if not xlsfn:
        warn('No XLSX parameter file found')
        return

    ap,sim,cam,Fwd = getParams(xlsfn,overrides,makeplot,progms,verbose=verbose)
    cam = definecamind(cam,sim.nCutPix)

    for a in ap:
        #TODO assumes all are same distance apart
        x0true = (ap[a].loc['X0km',:][:-1] + 0.5*diff(ap[a].loc['X0km',:]))[tInd]
        E0true = (ap[a].loc['E0',:][:-1]   + 0.5*diff(ap[a].loc['E0',:]))[tInd]

        analyseres(sim,cam,
                   x, xp, Phifwd, Phidict, drn, dhat,
                   vlim, x0true,E0true,makeplot,
                   progms, verbose=verbose)

#%%
    for ti,t in enumerate(tInd):
        if 'fwd' in makeplot:
            plotB(drn[ti],sim.realdata,cam,vlim['b'],t,2893,makeplot,'$b_{fwd',progms,verbose)

            plotJ(sim,Phifwd[...,ti],x,xp,Phidict[ti]['EK'],None,
                  vlim['j'],vlim['p'][:2],t,makeplot,'phifwd',
                  '$\phi_{top}$ fwd diff. number flux',
                  None,None,progms,verbose)

            plotVER(sim,Pfwd[ti],x,xp,z,zp,vlim['p'],t,makeplot,'pfwd',
                    '$P$ fwd volume emission rate',
                    None,None,progms,verbose)

        if 'optim' in makeplot:
           # plotB(dhat[ti],sim.realdata,cam,vlim['b'],t,2894,makeplot,'$b_{est',progms,verbose)
            plotBcompare(sim,drn[ti],dhat[ti],cam,'best',None,vlim['b'],t,2894,makeplot,progms,verbose)

            plotJ(sim,Phidict[ti]['x'],x,xp,Phidict[ti]['EK'],None,
                  vlim['j'],vlim['p'][:2],t,makeplot,'phiest',
                  '$\hat{\phi}_{top}$ estimated diff. number flux',
                  None,None,progms,verbose)

            plotVER(sim,Pest[ti],x,xp,z,zp,vlim['p'],t,makeplot,'pest',
                    '$\hat{P}$ estimated volume emission rate',
                    None,None,progms,verbose)

def findxlsh5(h5path):
    h5path = expanduser(h5path)

    if isfile(h5path):
        h5list = [h5path]
        xlsfn = glob(join(dirname(h5path),'*.xlsx'))
    else:
        h5list = glob(join(h5path,'dump_*.h5'))
        h5list.sort()
        xlsfn = glob(join(h5path,'*.xlsx'))

    if xlsfn: xlsfn = xlsfn[0]

    return h5list,xlsfn

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="load HiST output and plot/analyse")
    p.add_argument('h5path',help='path containing dump.h5 outputs (Hist output)')
    p.add_argument('-m','--makeplot',help='plots to make',default=[],nargs='+')
    p = p.parse_args()

    h5list,xlsfn = findxlsh5(p.h5path)

    readresults(h5list,xlsfn,None,p.makeplot)

    show()