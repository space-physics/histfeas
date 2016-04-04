#!/usr/bin/env python3
"""
To generate inputs for this program, run main_hist.py with
-m h5
option
"""
from pathlib import Path
import logging
import h5py
from numpy import asarray,diff
from matplotlib.pyplot import show
#
from histutils.findnearest import find_nearest
from .analysehst import analyseres
from .sanityCheck import getParams
from .plotsnew import plotoptim,plotfwd
from .observeVolume import definecamind

def readresults(h5list,xlsfn,vlim,x1d,overrides,makeplot,verbose=0):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]; Pest=[]; Pfwd=[]; ut1_unix=[]
#%%
    if not h5list:
        raise ValueError('no HDF5 files found')

    for h5 in h5list:
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

#                d = f['/best/bfit'].value
#                #TODO temp hack
#                d[512:] = d[512:][::-1]
#                dhat.append(d)
#
#                d = f['/best/braw'].value
#                d[512:] = d[512:][::-1]
#                drn.append(d)

                dhat.append(f['/best/bfit'].value)
                drn.append(f['/best/braw'].value)

                angle_deg = f['/best/angle'].value #NOTE: by definition, same angles for all time steps-the camera is not moving!

                try: #realdata
                    ut1_unix.append(f['/best/ut1_unix'].value)
                except KeyError: #simultation, not real data
                    pass

            except KeyError as e:
                raise KeyError('It seems that data inversion did not complete? Or at least it was not written  {}'.format(e))

#%%

    odir = h5.parent / 'reader'
    odir.mkdir(parents=True,exist_ok=True)

    try:
        Phifwd = asarray(Phifwd).transpose(1,2,0)
    except ValueError: #realdata
        pass

    if not xlsfn:
        raise ValueError('No XLSX parameter file found')

    ap,sim,cam,Fwd = getParams(xlsfn,overrides,makeplot,odir)
    cam = definecamind(cam,sim.nCutPix)
#%% load original angles of camera
    ut1_unix = asarray(ut1_unix)
    for i,C in enumerate(cam):
        C.angle_deg = angle_deg[i,:]
        C.tKeo = ut1_unix[:,i]
#%% load args if they exist
    for a in ap:
        #TODO assumes all are same distance apart
        x0true = (ap[a].loc['X0km',:][:-1] + 0.5*diff(ap[a].loc['X0km',:]))#[tind]
        E0true = (ap[a].loc['E0',:][:-1]   + 0.5*diff(ap[a].loc['E0',:]))  #[tInd]

        analyseres(sim,cam,
                   x, xp, Phifwd, Phidict, drn, dhat,
                   vlim, x0true,E0true,makeplot, odir)

#%% plots
    for i in range(len(drn)):
        try:
            Jxi = find_nearest(x,x1d[i])[0]
        except TypeError:
            Jxi = None

        try: #simulation
            pf = Pfwd[i]
            phif = Phifwd[...,i]
        except IndexError: #realdata
            pf = None;  phif = None

        if 'fwd' in makeplot:
            plotfwd(sim,cam,drn[i],x,xp,z,zp,
                    Pfwd[i],Phifwd[...,i],Phidict[i],Jxi,vlim,i,makeplot,odir,
                    doSubplots=True)

        if 'optim' in makeplot:
            plotoptim(sim,cam,drn[i],dhat[i],'best',pf,phif,Jxi,
                      Pest[i],Phidict[i],x,xp,z,zp,vlim,i,makeplot,odir,
                      doSubplots=True)


def findxlsh5(h5path):
    h5path = Path(h5path).expanduser()

    if h5path.is_file():
        h5list = [h5path]
        xlsfn = sorted(h5path.parent.glob('*.xls*'))
    elif h5path.is_dir():
        h5list = sorted(h5path.glob('dump*.h5'))
        xlsfn =  sorted(h5path.glob('*.xls*'))
    else:
        raise ValueError('no path or file at {}'.format(h5path))

    if isinstance(xlsfn,list):
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