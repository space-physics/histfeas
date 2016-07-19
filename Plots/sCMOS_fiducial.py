#!/usr/bin/env python
"""
creates sCMOS fiducials for TGARS 2015 Reconstruction of fine scale auroral Dyanmics paper figure
The AVI reader accounts for the full annotated Matlab axes frame and indexes accordingly
"""
from numpy import arange,datetime64,timedelta64
import cv2
from matplotlib.pyplot import show
#
from cvutils.getaviprop import getaviprop
from histfeas import Path
from histfeas.fiducial import fiducial
#%% sCMOS Figure
path = '~/data/2011-03-01/optical'
vidfn= 'CMOS_110301_1006.avi'
cadence_sec = 0.02 # by inspection of video title
t0 = datetime64('2011-03-01T10:06:09.06')
dt = timedelta64(int(cadence_sec*1e3),'ms')

tlim = (datetime64('2011-03-01T10:06:06'),
        datetime64('2011-03-01T10:06:22'))
#%%
fn = Path(path).expanduser() / vidfn
finf = getaviprop(fn)

t = arange(t0, finf['nframe']*dt, dt)

tlist = t[(t>tlim[0]) & (t<tlim[1])]

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

vid = cv2.VideoCapture(str(fn))

for i,T in enumerate(tlist):
    vid.set(cv2.CAP_PROP_POS_FRAMES, i) # advance to requested frame number 0-based
    iact = int(vid.get(cv2.CAP_PROP_POS_FRAMES))
    assert iact == i,'got frame {} but requested frame {}'.format(iact,i)

    ret,img = vid.read()
#%% time label -- elapsed time since first frame assuming uniform frame rate
    tstr = str(T-tlist[0])[:-12] + 'ms'
#%% plot!
    outfn = 'anno_{}.png'.format(T)
    fiducial(img, xycrop[0], xycrop[1], outfn, rings, rays,
             T,tstr, magzenith_xy, wh0, ringmult, axlim)

    rings=rays=False

show()