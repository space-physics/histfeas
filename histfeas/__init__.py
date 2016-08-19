from tempfile import mkdtemp
from sys import argv

try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    from pathlib2 import Path
#%%

def hist_figure(P):
    from histfeas.main_hist import doSim

    P['overrides']['rootdir'] = Path(__file__).parents[1]

    P['makeplot'] =['fwd','optim','png','h5']

    print('running HiSTfeas program -- will write png and h5 to {}'.format(P['outdir']))

    doSim(P)


def userinput(ini=None):
    from argparse import ArgumentParser
    p = ArgumentParser(description='flaming figure plotter')
    p.add_argument('ini',nargs='?',default=ini)
    p.add_argument('--load',help='load without recomputing',action='store_true')
    p.add_argument('-m','--makeplot',help='plots to make',default=[],nargs='+')
    p.add_argument('-L','--ell',help='compute projection matrix',action='store_true')
    p.add_argument('-v','--verbose',help='verbosity',action='count',default=0)
    p.add_argument('-f','--frames',help='time steps to use',nargs='+',type=int)
    p.add_argument('-o','--outdir',help='output directory',default = mkdtemp(dir='out'))
    p = p.parse_args()

    if not p.ini:
        raise RuntimeError('you must specify an .ini file')

#%% must occur first
    import matplotlib
    matplotlib.use('Agg')
    print(matplotlib.get_backend())
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
     'cmd': ' '.join(argv),
    }
#%%
    if p.frames is None or len(p.frames) not in (2,3):
        itime = p.frames
    elif len(p.frames)==2:
        itime = range(p.frames[0],p.frames[1])
    elif len(p.frames)==3:
        itime = range(p.frames[0],p.frames[1],p.frames[2])

    P['timeinds'] = itime


    return P
