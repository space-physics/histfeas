#!/usr/bin/env python
"""
creates sCMOS fiducials for TGARS 2015 Reconstruction of fine scale auroral Dyanmics paper figure
The AVI reader accounts for the full annotated Matlab axes frame and indexes accordingly
"""
from numpy import arange
import cv2
from matplotlib.pyplot import show
#
from histfeas import Path
from histfeas.fiducial import fiducial
#%% sCMOS Figure
path = '~/data/2011-03-01/optical'
vidfn= 'CMOS_110301_1006.avi'
flist = range(100,102) # frame indices

cadence_sec = 0.02 # by inspection of video title

xycrop=(0,0) # normally (0,0) for uncropped original video
#%% zero-based indexing from upper left of image, magnetic zenith
magzenith_xy = (1328,876)
#%% half-width, half-height of the first oval, that is, 5 degreees elevation = wh0[1]
geozenith_xy = (1242,297) #geographic zenith
wh0 = (1494-geozenith_xy[0], 531-geozenith_xy[1])

axlim = (magzenith_xy[0]-500, magzenith_xy[0]+500,
         magzenith_xy[1]+275, magzenith_xy[1]-275) #optional pixel indices to show (xmin,xmax, ymax,ymin)
ringmult = arange(2,12,2,dtype=float) # magnetic zenith angle to label with rings
# for first image only
rings=True
rays=False

path = Path(path).expanduser()

vid = cv2.VideoCapture(str(path/vidfn))

for f in flist:
    vid.set(cv2.CAP_PROP_POS_FRAMES,f) # advance to requested frame number 0-based
    fact = int(vid.get(cv2.CAP_PROP_POS_FRAMES))
    assert fact == f,'got frame {} but requested frame {}'.format(fact,f)

    ret,img = vid.read()
#%% time label -- elapsed time since first frame assuming uniform frame rate
    tstr = '{:.0f} ms'.format( (f-flist[0]) * cadence_sec * 1e3)
#%% plot!
    outfn = 'anno_{}.png'.format(f)
    fiducial(img, xycrop[0], xycrop[1], outfn, rings, rays,
             tstr, magzenith_xy, wh0, ringmult, axlim)

    rings=rays=False

show()