#!/usr/bin/env python3
"""
generates circular markings of az/el on videos
 montage anno_flame10061*.png -trim -tile 4x1 -geometry +1+0  out.png
"""
from __future__ import division
from six import PY2
from . import Path
from matplotlib.patches import Ellipse
from matplotlib.pyplot import figure, Axes,close
import matplotlib
matplotlib.rcParams.update({'font.family':'sans-serif',
                           'font.sans-serif':'Arial',
                           'text.usetex':True})
from scipy.misc import imread
from scipy.io import loadmat
from numpy import array,cos,sin,radians,ndarray
#
from GeoData.plotting import plotazelscale

RINGCALDEG = 5.
DPI = 75

if PY2: FileNotFoundError = OSError

def ccdfid(imgfn,calfn):
    imgfn = Path(imgfn).expanduser()
    calfn = Path(calfn).expanduser()
#%% load cal data
    cal = loadmat(str(calfn))
#%%
    opt = imread(str(imgfn))
    plotazelscale(opt,cal['az'],cal['el'])

def fiducial(img,xcrop,ycrop,outfn,rings,rays,t,pstr,oxyfull,wh0, ringmult,axlim=None):

    """
    img: geographic zenith, full frame video frame including axes/margins
    """
    assert isinstance(img,ndarray) and img.ndim in (2,3),'I expect a single image, greyscale or RGB'
    assert isinstance(xcrop,int)
    assert isinstance(ycrop,int)
    assert isinstance(rings,bool)
    assert isinstance(rays,bool)
#%% shift to account for post-processing crop of video (by me)
    oxy = [oxyfull[0]-xcrop, oxyfull[1]-ycrop]
#%% init image
    fg = figure()
    ax = Axes(fg, [0., 0., 1., 1.]) #square
    ax.set_axis_off()
    fg.add_axes(ax)
    if axlim:
        ax.axis(axlim)
    ax.imshow(img,cmap='gray',origin='upper')
    #ax.autoscale(True) #scale to image, not rings
    #%% plot rings
    # algorithm assumes uniform spaced rings
    if rings:
        xyr1 = array((wh0[0]*ringmult[0]/RINGCALDEG-xcrop,
                      wh0[1]*ringmult[0]/RINGCALDEG-ycrop)) #by inspection
        ringincr = wh0[1]*ringmult[0]/RINGCALDEG #every 2 degrees
        for m in ringmult:
            xyr = m*xyr1 + 1.

            ax.add_patch(Ellipse(xy=oxy,
                         width=xyr[0], height=xyr[1],
                         fc='none',ec='white',lw=1))

            angtxt(ringincr*m/2., oxy,'{:.0f}$^\circ$'.format(m), ax, wh0)
    #%% plot rays
    magzenel = 77.5
    magzenaz = 207.5
    # 1328, 876
    if rays:
        xyl = []
        #angdeg = [40,70,85.5,98.75,117]
        az180 = 109.65 #offset in degrees to 180 degree ray
        angdeg = az180 - (array([magzenaz]) - 180 )

        elr = wh0[1]/5*(90-magzenel)

        for a in radians(angdeg):
            #by equation of general ellipse
            xyl.append([wh0[0]/wh0[1]*elr*cos(a)+oxy[0], elr*sin(a)+oxy[1]])

        for r in xyl:
            ax.plot([r[0],oxy[0]], [r[1],oxy[1]],
                    color='red',linestyle='-',lw=1)
#%% annotate times
    ax.text(0.035,0.965,pstr,color='white',fontsize=50,
            va='top',ha='left', transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))

    ax.set_title(str(t), fontsize='x-large')
    #%% finalize image
    #ax.axis('off')
    #fg.tight_layout

    if outfn:
        print('writing {}'.format(outfn))
        fg.savefig(str(outfn),bbox_inches='tight',dpi=DPI)
        close(fg)


def angtxt(radput,oxy,txt,ax,wh0):
    angput = radians(30)
    ax.text(wh0[0]/wh0[1]*radput*cos(angput)+oxy[0],
            radput*sin(angput)+oxy[1],
            txt,color='white',fontsize=30,
            va='top',ha='left',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))
