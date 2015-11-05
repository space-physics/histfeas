#!/usr/bin/env python3
"""
To generate inputs for this program, run main_hist.py with
-m h5
option
"""

from __future__ import division,absolute_import
from pathlib2 import Path
import logging
import re
import h5py
from os import makedirs
from numpy import asarray,diff
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='poster', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'})
#
from histutils.findnearest import find_nearest
from .analysehst import analyseres
from .sanityCheck import getParams
from .plotsnew import plotoptim,plotfwd
from .observeVolume import definecamind

def readresults(h5list,xlsfn,vlim,x1d,overrides,makeplot,verbose=0):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]; Pest=[]; Pfwd=[]
    tInd = [];
#%%
    nt = len(h5list)
    if nt==0:
        logging.warning('no HDF5 files found, ending.')
        return

    for h5 in h5list:
        tInd.append(re.compile(r'(\d+)$').search(h5.stem).group(1)) #NOTE assumes trailing time integers

        logging.debug('tind {}  reading {}'.format(tInd[-1],h5))

        assert h5.is_file()
        with h5py.File(str(h5),'r',libver='latest') as f:
            try: #simulation
                Phifwd.append(f['/phifwd/phi'].value)
                Pfwd.append(f['/pfwd/p'].value)
            except KeyError: #real data
                pass

            try:
                Phidict.append({'x':f['/phiest/phi'].value,
                                'EK':f['/phiest/Ek'].value,
                                'EKpcolor':f['/phiest/EKpcolor'].value})

                Pest.append(f['/pest/p'].value)

                x  = f['/pest/x'].value #same for all in directory
                xp = f['/pest/xp'].value
                z = f['/pest/z'].value
                zp = f['/pest/zp'].value

                dhat.append(f['/best/bfit'].value)
                drn.append(f['/best/braw'].value)
            except KeyError as e:
                raise KeyError('It seems that data inversion did not complete? Or at least it was not written  {}'.format(e))
#%%

    progms = h5.parent / 'reader'

    makedirs(progms,exist_ok=True)


    try:
        Phifwd = asarray(Phifwd).transpose(1,2,0)
    except ValueError: #realdata
        pass

    if not xlsfn:
        logging.error('No XLSX parameter file found')
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
        try:
            Jxi = find_nearest(x,x1d[ti])[0]
        except:
            Jxi = None

        try:
            pf = Pfwd[ti]
            phif = Phifwd[...,ti]
        except:
            pf = None
            phif = None

        if 'fwd' in makeplot:
            plotfwd(sim,cam,drn[ti],x,xp,z,zp,
                    Pfwd[ti],Phifwd[...,ti],Phidict[ti],Jxi,vlim,t,makeplot,None,progms,verbose)

        if 'optim' in makeplot:
            plotoptim(sim,cam,drn[ti],dhat[ti],'best',pf,phif,Jxi,
                      Pest[ti],Phidict[ti],x,xp,z,zp,vlim,t,makeplot,None,progms,verbose)


def findxlsh5(h5path):
    h5path = Path(h5path).expanduser()

    if h5path.is_file():
        h5list = [h5path]
        xlsfn = sorted(h5path.parent.glob('*.xls*'))
    else:
        h5list = sorted(h5path.glob('dump_*.h5'))
        xlsfn =  sorted(h5path.glob('*.xls*'))

    if xlsfn:
        xlsfn = xlsfn[0]

    return h5list,xlsfn

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description="load HiST output and plot/analyse")
    p.add_argument('h5path',help='path containing dump.h5 outputs (Hist output)')
    p.add_argument('-m','--makeplot',help='plots to make',default=[],nargs='+')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p = p.parse_args()

    h5list,xlsfn = findxlsh5(p.h5path)

    readresults(h5list,xlsfn,vlim=None,Jxi=None,overrides=None,
                makeplot=p.makeplot,verbose=p.verbose)

    show()