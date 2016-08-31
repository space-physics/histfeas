#!/usr/bin/env python
"""
creates EMCCD fiducials for TGARS 2015 Reconstruction of fine scale auroral Dyanmics paper figure
"""
from numpy import array,rot90
from scipy.ndimage import imread
from matplotlib.pyplot import show
#
from histfeas import Path
from histfeas.fiducial import fiducial
#%% EMCCD Figure
path = '~/data/2007-03/optical'
ccdflist = ('1132.png','1147.png','1162.png','1177.png')
ccdcal= '~/data/CMOS/X1387_03_23_2007_031836.mat'

xycrop=(0,0)
#    ccdfid(path/ccdflist[0],ccdcal)

oxyfull = [222,190]
wh0 = (140,140)
pstr=('0s','0.5s','1.0s','1.5s')
lblring = array((89,87,85))
ringmult=90-lblring

rings=(True,False,False,False)
rays=(False,False,False,False) # set first True, and oxyfull to geographic zenith to plot ray pointing to magnetic zentih

path = Path(path).expanduser()

for f,ring,ray,p in zip(ccdflist,rings,rays,pstr):
    imgfn = path/f
    outfn = path/('anno_' + f)

    try:
        img = imread(str(imgfn))
        img = rot90(img,-1)
    except FileNotFoundError:
        continue

    fiducial(img, xycrop[0], xycrop[1], outfn, ring, ray, p,
             oxyfull,wh0,lblring,ringmult,wh0)

show()