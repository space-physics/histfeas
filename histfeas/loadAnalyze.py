#!/usr/bin/env python
"""
To generate inputs for this program, run main_hist.py with
-m h5
option
"""
from __future__ import division
from . import Path
import h5py
from numpy import asarray,diff
from matplotlib.pyplot import show,close
#
from .analysehst import analyseres
from . import getParams
from .plotsnew import plotoptim,plotfwd
from .observeVolume import definecamind

def readresults(h5list, P):
    Phifwd =[]; Phidict =[]; dhat=[]; drn=[]; Pest=[]; Pfwd=[]; ut1_unix=[]
#%%
    if not h5list:
        print('no HDF5 files found from your analysis run.')
        return

    if len(h5list)>500:
        print('loading {} files from {}'.format(len(h5list),h5list[0].parent))

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
                                'EKpcolor':f['/phiest/EKpcolor'].value,
                                'gx0':f['/phiest/gx0'].value,
                                'gE0':f['/phiest/gE0'].value,})

                Pest.append(f['/pest/p'].value)

                x  = f['/pest/x'].value #same for all in directory
                xp = f['/pest/xp'].value
                z = f['/pest/z'].value
                zp = f['/pest/zp'].value

                dhat.append(f['/best/bfit'].value)
                drn.append(f['/best/braw'].value)

                angle_deg = f['/best/angle'].value #NOTE: by definition, same angles for all time steps-the camera is not moving!

                try: #realdata
                    ut1_unix.append(f['/best/ut1_unix'].value)
                except KeyError: #simulation, not real data
                    pass

            except KeyError as e:
                raise KeyError('It seems that data inversion did not complete? Or at least it was not written  {}'.format(e))
#%%

    odir = h5.parent / 'reader'
    print('writing output to {}'.format(odir))
    odir.mkdir(parents=True,exist_ok=True)

    if Phifwd: #sim data
        Phifwd = asarray(Phifwd).transpose(1,2,0) #result: Nenergy x Nx x Ntime

    arc,sim,cam,Fwd = getParams(P)
    cam = definecamind(cam)
#%% load original angles of camera
    ut1_unix = asarray(ut1_unix)
    for i,C in enumerate(cam):
        if C.usecam:
            C.angle_deg = angle_deg[i,:]
            C.tKeo = ut1_unix[:,i]
#%% load args if they exist
    for k,a in arc.items():
        """
        TODO: assumes all are same distance apart (zero accel)
        NOTE: Assumption is that arc is moving smoothly with zero acceleration,
        so that mean position for each time step is = x[:-1] + 0.5*diff(x) as below.
        """
        if a.zshape in ('flat','impulse'):
            x0true = a.X0km[:-1]
            E0true = a.E0[:-1]
        else:
            x0true = (a.X0km[:-1] + 0.5*diff(a.X0km))
            E0true = (a.E0[:-1]   + 0.5*diff(a.E0))

        analyseres(sim,cam, x, xp, Phifwd, Phidict, drn, dhat,P, x0true,E0true)

#%% plots
    for i in range(len(drn)): #for each time, do...
        try: #simulation
            pf = Pfwd[i]
            phif = Phifwd[...,i]
        except IndexError: #realdata
            pf = None;  phif = None

        if 'fwd' in P['makeplot']:
            plotfwd(sim,cam,drn[i],x,xp,z,zp, pf,phif,Phidict[i],i,P,doSubplots=True)

        if 'optim' in P['makeplot']:
            plotoptim(sim,cam,drn[i],dhat[i],'best',pf,phif, Pest[i],Phidict[i],x,xp,z,zp,i,P,doSubplots=True)

        if 'show' in P['makeplot']:
            show()
        else:
            close('all')

    return Phifwd, Phidict

def findxlsh5(h5path):
    h5path = Path(h5path).expanduser()

    if h5path.is_file():
        h5list = [h5path]
        inifn = sorted(h5path.parent.glob('*.ini'))
    elif h5path.is_dir():
        h5list = sorted(h5path.glob('dump*.h5'))
        inifn =  sorted(h5path.glob('*.ini'))
    else:
        raise FileNotFoundError('no path or file at {}'.format(h5path))

    if not inifn:
        raise FileNotFoundError('no simulation ini found in {}'.format(h5path))

    inifn = inifn[0]



    return h5list,inifn

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
