#!/usr/bin/env python3
"""
generates circular markings of az/el on videos
 montage anno_flame10061*.png -trim -tile 4x1 -geometry +1+0  out.png
"""
from __future__ import division
from six import PY2
from . import Path
from matplotlib.patches import Ellipse
from matplotlib.pyplot import figure, Axes
import matplotlib
matplotlib.rcParams.update({'font.family':'sans-serif',
                           'font.sans-serif':'Arial',
                           'text.usetex':True})
from scipy.misc import imread
from scipy.io import loadmat
from numpy import array,cos,sin,radians,ndarray
#
from GeoData.plotting import plotazelscale

if PY2: FileNotFoundError = OSError

def ccdfid(imgfn,calfn):
    imgfn = Path(imgfn).expanduser()
    calfn = Path(calfn).expanduser()
#%% load cal data
    cal = loadmat(str(calfn))
#%%
    opt = imread(str(imgfn))
    plotazelscale(opt,cal['az'],cal['el'])

def fiducial(img,xcrop,ycrop,outfn,rings,rays,pstr,oxyfull,wh0,
             lblring,ringmult,axlim=None):

    """
    img: geographic zenith, full frame video frame including axes/margins
    """
    assert isinstance(img,ndarray) and img.ndim in (2,3),'I expect a single image, greyscale or RGB'
    assert isinstance(xcrop,int)
    assert isinstance(ycrop,int)
    assert isinstance(rings,bool)
    assert isinstance(rays,bool)
#%% setup rings
    xyr = []
    #zero-based indexing, geographic zenith
    #oxyfull = [1242,297]
#%% shift to account for post-processing crop of video (by me)
    oxy = [oxyfull[0]-xcrop, oxyfull[1]-ycrop]

    xyr.append(array([wh0[0]*2/5-xcrop, wh0[1]*2/5-ycrop])) #by inspection

#%% width,height of rest of ovals
    for i in ringmult:
        xyr.append(i*xyr[0]+1) #multiple of radius
#%% init image
    fg = figure()
    ax = Axes(fg, [0., 0., 1., 1.])
    ax. set_axis_off()
    fg.add_axes(ax)
    if axlim:
        ax.axis(axlim)
    ax.imshow(img,cmap='gray',origin='upper')
    #ax.autoscale(True) #scale to image, not rings
    #%% plot rings
    if rings:
        ringincr = wh0[1]*2/5. #every 2 degrees
        for c in xyr:
            ax.add_patch(Ellipse(xy=oxy,
                             width=2*c[0], height=2*c[1],
                             fc='none',ec='white',lw=1))
        #%% annotate rings
        for m,r in zip(ringmult,lblring):
            angtxt(ringincr*m, oxy,'{}$^\circ$'.format(r), ax, wh0)

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
    #%% finalize image
    #ax.axis('off')
    #fg.tight_layout

    #ax.set_ylim((500,1000))
    if outfn:
        print('writing {}'.format(outfn))
        fg.savefig(str(outfn),bbox_inches='tight',dpi=300)


def angtxt(radput,oxy,txt,ax,wh0):
    angput = radians(30)
    ax.text(wh0[0]/wh0[1]*radput*cos(angput)+oxy[0],
            radput*sin(angput)+oxy[1],
            txt,color='white',fontsize=30,
            va='top',ha='left',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))
