from pathlib import Path
import subprocess,os
import stat
from tempfile import mkdtemp
from argparse import ArgumentParser
from sys import argv
import logging
from configparser import ConfigParser
from shutil import copy2
from pymap3d.vincenty import vdist
#%%
from histutils.camclass import Cam
from histutils import splitconf
from .simclass import Sim
from .arcclass import Arc,getntimes
#%%

def hist_figure(P):
    from .main_hist import doSim # KEEP in this function to avoid ImportError
    print(f'running HiSTfeas program -- will write png and h5 to {P["outdir"]}')
    doSim(P)


def userinput(ini:Path=None, outdir:Path=None) -> dict:

    p = ArgumentParser(description='flaming figure plotter')
    p.add_argument('ini',help='.ini config file',nargs='?',default=ini)
    p.add_argument('outdir',help='output directory',nargs='?',default = outdir)

    p.add_argument('-g','--fwdguess',help='feed minimizer fwd answer. true | randn stddev |',nargs='+')
    p.add_argument('--cx',help='override cam positions (must specify all)',nargs='+',type=float)
    p.add_argument('--iter',help='number of data inversion iterations',type=int)
    p.add_argument('--fitm',help='fit method')

    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-m','--makeplot',help='plots to make',default=['realvid','fwd','optim'],nargs='+')
    p.add_argument('-L','--ell',help='compute projection matrix',action='store_true')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',nargs='+',type=int)
    p = p.parse_args()

    if not p.ini:
        raise FileNotFoundError('you must specify an .ini file')
#%%
    P = {'ini':Path(p.ini).expanduser(),
     'load':p.load,
     'makeplot':p.makeplot,
     'ell':p.ell,
     'verbose':p.verbose,
     'overrides': {},
     'cmd': ' '.join(argv),
     'gitrev': subprocess.check_output(['git','rev-parse','--short','HEAD']).decode('utf8').strip('\n')
    }

#%% directory handling
    if not p.load and not p.outdir:
        p.outdir = mkdtemp(prefix=P['ini'].stem,dir='out')

    if P['load'] and not p.outdir:
        if P['ini'].is_file():
            p.outdir = P['ini'].parent
        elif P['ini'].is_dir():
            p.outdir = P['ini']
        else:
            raise FileNotFoundError('nothing found at {}'.format(P['ini']))

    P['outdir'] = Path(p.outdir).expanduser()
    P['outdir'].mkdir(parents=True,exist_ok=True)
#%%
    P['overrides']['rootdir'] = Path(__file__).parents[1]

    P['overrides']['fwdguess'] = p.fwdguess
    P['overrides']['camx'] = p.cx
    P['overrides']['fitm'] = p.fitm
    P['overrides']['niter'] = p.iter
#%%
    if p.frames is None or len(p.frames) not in (2,3):
        itime = p.frames
    elif len(p.frames)==2:
        itime = range(p.frames[0],p.frames[1])
    elif len(p.frames)==3:
        itime = range(p.frames[0],p.frames[1],p.frames[2])

    P['timeinds'] = itime

    return P


def getParams(P:dict):
#%% first copy .ini file readonly to output dir for future reference
    if P['outdir'] is not None and not P['load']:  #running new simulation
        try: # so we can overwrite existing file when user wants to reuse output directory
            os.chmod(P['outdir']/P['ini'].name, stat.S_IROTH+stat.S_IRUSR+stat.S_IRGRP+stat.S_IWUSR)
        except FileNotFoundError:
            pass

        copy2(P['ini'], P['outdir'])
      # make the copied-to-output-directory .ini file readonly to help avoid mistakes when loading and replotting
        os.chmod(P['outdir'] / P['ini'].name, stat.S_IROTH+stat.S_IRUSR+stat.S_IRGRP)
#%% read .ini
    xl = ConfigParser(allow_no_value=True, empty_lines_in_values=False,
                      inline_comment_prefixes=(';'), strict=True)

    xl.read(P['ini'])
# %% read plot parameters from ini
    if 'plot' in xl:
        P = plotstuffer(xl['plot'],P)
    else:
        P = plotstuffer(None,P)
#%% read arcs (if any)
    arc = setupArc(xl)
    logging.info(f'# of observer time steps in {P["ini"]}: {arc["ntimeslice"]}')
#%% class with parameters and function
    sim = Sim(xl,arc,P)
#%% grid setup
    Fwd = sim.setupFwdXZ(xl)
#%% setup cameras
    cam = setupCam(sim, xl['cam'], Fwd['z'][-1], P)

    if cam[0].x_km is None: #defined simulated camera location in lat/lon
        cam = cam0dist(cam)

    #store x,z in sim
    ellname=sim.getEllHash(xl, [C.x_km for C in cam if C.usecam],
                               [C.alt_m/1000. for C in cam if C.usecam])
    #will try to load this and compute if needed. Will be copied to output directory too.
    sim.FwdLfn = sim.rootdir/'precompute' / ellname

    # make the simulation time step match that of the fastest camera
    sim.kineticsec = min([C.kineticsec for C in cam if C.usecam])


    logging.info('fwd model voxels:\n'
          f'B_perp: N={Fwd["sx"]}   B_parallel: M={Fwd["sz"]}')
#%% init variables
    return arc,sim,cam,Fwd,P


def cam0dist(cam):
    """
    find the x-coordinate (along B-perp) of each camera (can't do this inside camclass.py)
    """
    cam[0].x_km = 0. # NOTE arbitrarily picks the first camera x=0km
    for C in cam[1:]:
        if C.usecam:
            C.x_km = vdist(cam[0].lat,cam[0].lon,C.lat,C.lon)[0] / 1e3

    return cam


def setupArc(xl) -> dict:
    """
    in the .ini file, sections [arc] using csv for each arc
    will create synthetic auroral arcs
    """
    arc = {}
    ntimeslice=None

    for s in xl:
        if not s.startswith('arc'):
            continue

        logging.debug(f'configuring {s}')
        texp = getntimes(xl[s]['texp']) # last time is blended with 2nd to last time
        if ntimeslice is not None:
            assert len(texp) == ntimeslice+1, 'for now, all Arcs must have same number of times (columns)'

        # TODO assert all arcs have same time length
        ntimeslice = texp.size-1  # MUST be -1 to allow for range() compat and interp1() compat

        arc[s] = Arc(xl[s], texp)

    arc['ntimeslice'] = ntimeslice

    return arc


def setupCam(sim, cp, zmax:float, P:dict):
    cam = []

    if sim.camxreq is not None:
        logging.warning(f'overriding camera x-loc with {sim.camxreq}')
        for C,cx in zip(sim.camnames,sim.camxreq): #enumerate as in general camera 0 may not be used
            cam.append(Cam(sim, cp, C, zmax, xreq=cx,
                           makeplot=P['makeplot'], verbose=P['verbose']))
    else:
        for C in sim.camnames:
            cam.append(Cam(sim, cp, C, zmax,
                           makeplot=P['makeplot'], verbose=P['verbose']))

    assert len(cam)>0,'0 cameras are configured, Nothing to do.'

    return cam


def plotstuffer(sp, P:dict) -> dict:
    """
    these have no impact on simulation calculations, they are just plotting bounds
    """
    if sp is not None and 'x1d' in sp:
        P['x1d'] = splitconf(sp,'x1d')

    if 'vlim' not in P:
        P['vlim'] = {}

    for p in ('p','p1d','j','j1d','b','x','z'):
        P['vlim'][p] = splitconf(sp,p,fallback=(None,None))

    #print(P['vlim'])
    return P
