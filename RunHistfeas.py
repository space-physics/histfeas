#!/usr/bin/env python3
"""
Example command line interfacer for HIST feasibility
Michael Hirsch
"""

import matplotlib as mpl
from matplotlib.pyplot import show
#
from histfeas import userinput
from histfeas.main_hist import doSim

print('matplotlib backend:',mpl.get_backend(),' version', mpl.__version__)

if __name__ == '__main__':
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    #p = ArgumentParser(description='analyzes HST data and makes simulations')
   # p.add_argument('infile',help='.xls filename with simulation parameters')
   # p.add_argument('outdir',help='directory for output')
   # p.add_argument('-v','--verbose',help='set debugging verbosity e.g. -v -vv -vvv',action='count',default=0)
   # p.add_argument('-m','--makeplot',help='list plots you want made',nargs='+',default=['']) #None gave type errors in generators & list comp.
  #  p.add_argument('-f','--frames',help='START STOP STEP of time indices',nargs=3,type=int)
  #  p.add_argument('--profile',help='Profile performance of program (development/debug only)',action='store_true')
  #  p.add_argument('-c','--cam',help='zero-indexed cameras to use (overrides XLS)',nargs='+',default=[None],type=int)
 #   p.add_argument('--cx',help='override cam positions (must specify all)',nargs='+',type=float)
 #   p.add_argument('--influx',help='flux .h5 file to use (overrides XLS)')
 #   p.add_argument('--x1d',help='required location [km] of x for 1-D flux plots',nargs='+',default=[None],type=float)
#    p.add_argument('--saveall',help='saves "all" variables to disk, useful for registration case tracking',action='store_true')
 #   p.add_argument('-a','--anim',help='animate plots (crudely)',type=float)
 #   p.add_argument('--filter',help='optical filter choices: bg3   none')
#    p.add_argument('--minev',help='minimum beam energy to include in fwd model',type=float)
#    p.add_argument('-g','--fwdguess',help='feed minimizer fwd answer. true | randn stddev |',nargs='+',default=[None])
#    p.add_argument('--fitm',help='override fit (minimization) method')
#    p.add_argument('--vlim',help='xmin xmax zmin zmax pmin pmax   limits for VER plots',type=float,nargs=6,default=[None]*6)
#    p.add_argument('--jlim',help='MIN MAX flux limits for diff num flux plots',type=float,nargs=2,default=(None,)*2)
#    p.add_argument('--blim',help='MIN MAX flux limits for brightness plots',type=float,nargs=2,default=(None,)*2)
#    p.add_argument('-L','--ell',help='force recomputation of sparse matrix L',action='store_true')
#    p = p.parse_args()

    P = userinput()
    doSim(P)

    show()
