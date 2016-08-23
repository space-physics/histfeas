import matplotlib
matplotlib.use('Agg')
print(matplotlib.get_backend())
#%%
try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    from pathlib2 import Path
#%%
from tempfile import mkdtemp
from argparse import ArgumentParser
from sys import argv
import logging
from six.moves.configparser import ConfigParser
from shutil import copy2
from geopy.distance import vincenty
#from numpy import arange, fromstring
#%%
from histutils.camclass import Cam
from .simclass import Sim
from .arcclass import Arc,getntimes
#%%
def hist_figure(P):
    from .main_hist import doSim # KEEP in this function to avoid ImportError
    P['makeplot'] =['fwd','optim','png','h5']

    print('running HiSTfeas program -- will write png and h5 to {}'.format(P['outdir']))

    doSim(P)


def userinput(ini=None,outdir=None):
    if ini:
        ini = Path(ini)

    if outdir is None:
        outdir = mkdtemp(prefix=ini.stem,dir='out')

    p = ArgumentParser(description='flaming figure plotter')
    p.add_argument('ini',help='.ini config file',nargs='?',default=ini)
    p.add_argument('outdir',help='output directory',nargs='?',default = outdir)
    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-m','--makeplot',help='plots to make',default=[],nargs='+')
    p.add_argument('-L','--ell',help='compute projection matrix',action='store_true')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',nargs='+',type=int)
    p = p.parse_args()

    if not p.ini:
        raise RuntimeError('you must specify an .ini file')
#%% now the other matplotlib imports
    import seaborn as sns
    sns.color_palette("cubehelix")
    sns.set(context='paper', style='whitegrid',font_scale=2,
            rc={'image.cmap': 'cubehelix_r'})
#%%
    P = {'ini':Path(p.ini).expanduser(),
     'load':p.load,
     'makeplot':p.makeplot,
     'ell':p.ell,
     'verbose':p.verbose,
     'outdir':Path(p.outdir).expanduser(),
     'overrides': {},
     'cmd': ' '.join(argv),
     'vlim':{'p':(None,None), 'p1d': (None,None),
             'j':(None,None), 'j1d': (None,None),
             'b':(None,None),
             'x':(None,None), 'z':(None,None)}
    }

    P['overrides']['rootdir'] = Path(__file__).parents[1]
#%%
    if p.frames is None or len(p.frames) not in (2,3):
        itime = p.frames
    elif len(p.frames)==2:
        itime = range(p.frames[0],p.frames[1])
    elif len(p.frames)==3:
        itime = range(p.frames[0],p.frames[1],p.frames[2])

    P['timeinds'] = itime

    return P



def getParams(P):
    if P['outdir'] is not None:
        copy2(str(P['ini']),str(P['outdir']))
#%% read .ini
    try:
        xl = ConfigParser(allow_no_value=True, inline_comment_prefixes=('#'), strict=True)
    except TypeError: #py27
        xl = ConfigParser(allow_no_value=True)
    xl.read(str(P['ini']))
#%% read arcs (if any)
    arc,ntimeslice = setupArc(xl)
    logging.info('# of observer time steps in {}: {}'.format(P['ini'],ntimeslice))
#%% class with parameters and function
    sim = Sim(xl,arc,ntimeslice,P)
#%% grid setup
    Fwd = sim.setupFwdXZ(xl)
#%% setup cameras
    cam = setupCam(sim,xl['cam'],Fwd['z'][-1],P['makeplot'])

    if sim.realdata:
        #find the x-coordinate (along B-perp) of each camera (can't do this inside camclass.py)
        cam[0].x_km = 0 # NOTE arbitrarily picks the first camera x=0km
        for C in cam[1:]:
            if C.usecam:
                C.x_km = vincenty((cam[0].lat,cam[0].lon),(C.lat,C.lon)).kilometers

    #store x,z in sim
    ellname=sim.getEllHash(xl, [C.x_km for C in cam],[C.alt_m/1000. for C in cam])
    #will try to load this and compute if needed. Will be copied to output directory too.
    sim.FwdLfn = sim.rootdir/'precompute' / ellname

    # make the simulation time step match that of the fastest camera
    sim.kineticsec = min([C.kineticsec for C in cam if C.usecam])


    logging.info('fwd model voxels:\n'
          'B_perp: N={}   B_parallel: M={}'.format(Fwd['sx'],Fwd['sz']))
#%% init variables
    return arc,sim,cam,Fwd

def setupArc(xl):
    """
    in the .ini file, sections [arc] using csv for each arc
    will create synthetic auroral arcs
    """
    arc = {}
    ntimeslice=None

    for s in xl.sections(): # for py27
        # TODO assert all arcs have same time length
        if s.startswith('arc'):
            if ntimeslice is not None and getntimes(xl[s]) != ntimeslice:
                raise ValueError('for now, all Arcs must have same number of times (columns)')
            texp = getntimes(xl[s]['texp']) # last time is blended with 2nd to last time

            arc[s] = Arc(xl[s], texp)

    return arc, texp.size-1 # MUST be -1 to allow for range() compat and interp1() compat

def setupCam(sim,cp,zmax,makeplot):
    cam = []

    if sim.camxreq[0] is not None:
        logging.warning('overriding camera x-loc with {}'.format(sim.camxreq))
        for C,cx in zip(sim.camnames,sim.camxreq): #enumerate as in general camera 0 may not be used
            cam.append(Cam(sim, cp, C, zmax,makeplot))
    else:
        for C in sim.camnames:
            cam.append(Cam(sim, cp, C, zmax,makeplot))

    assert len(cam)>0,'0 cams are configured, Nothing to do.'

    return cam