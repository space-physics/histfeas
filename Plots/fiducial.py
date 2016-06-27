#!/usr/bin/env python3
"""
generates circular markings of az/el on videos
 montage anno_flame10061*.png -trim -tile 4x1 -geometry +1+0  out.png
"""
from __future__ import division,absolute_import
from pathlib import Path
from matplotlib.patches import Ellipse
from matplotlib.pyplot import figure,show, Axes
import matplotlib
matplotlib.rcParams.update({'font.family':'sans-serif',
                           'font.sans-serif':'Arial',
                           'text.usetex':True})
from scipy.misc import imread
from scipy.io import loadmat
from numpy import array,cos,sin,radians,rot90
#
from GeoData.plotting import plotazelscale

def ccdfid(imgfn,calfn):
    imgfn = Path(imgfn).expanduser()
    calfn = Path(calfn).expanduser()
#%% load cal data
    cal = loadmat(str(calfn))
#%%
    opt = imread(str(imgfn))
    plotazelscale(opt,cal['az'],cal['el'])

def fiducial(imgfn,xcrop,ycrop,outfn,rings,rays,pstr,oxyfull,wh0,
             lblring,ringmult,axlim=None):

    # geographic zenith, full frame video including axes/margins
    #%% setup rings
    xyr = []
    #zero-based indexing, geographic zenith
    #oxyfull = [1242,297]

    #shift to account for post-processing crop of video (by me)
    oxy = [oxyfull[0]-xcrop, oxyfull[1]-ycrop]

    xyr.append(array([wh0[0]*2/5-xcrop, wh0[1]*2/5-ycrop])) #by inspection

    #width,height of rest of ovals
    for i in ringmult:
        xyr.append(i*xyr[0]+1) #multiple of radius
    #%% load image
    try:
        img = imread(str(imgfn))
        img = rot90(img,-1)
    except FileNotFoundError:
        return
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
            angtxt(ringincr*m, oxy,'{}$^\circ$'.format(r),ax)

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
    if outfn is not None:
        print('writing',outfn)
        fg.savefig(str(outfn),bbox_inches='tight',dpi=300)


def angtxt(radput,oxy,txt,ax):
    angput = radians(30)
    ax.text(wh0[0]/wh0[1]*radput*cos(angput)+oxy[0],
            radput*sin(angput)+oxy[1],
            txt,color='white',fontsize=30,
            va='top',ha='left',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='program to annotate cmos video avi')
    p.add_argument('--path',help='path of images',default='')
    p.add_argument('-c','--xycrop',help='origin of crop in x,y',nargs=2,default=(0,0),type=int)
    p.add_argument('-f','--file',help='name of single file to process')
    p.add_argument('-m','--magzen',help='azimuth from north (deg.) of magnetic zenith in image coordinate system',type=float)
    a = p.parse_args()

    path = Path(a.path).expanduser()

    rings=(True,False,False,False)
    rays=(False,False,False,False) # set first True, and oxyfull to geographic zenith to plot ray pointing to magnetic zentih

#%% EMCCD Figure 1
    ccdflist = ('1132.png','1147.png','1162.png','1177.png')
    ccdcal= '~/data/CMOS/X1387_03_23_2007_031836.mat'
#    ccdfid(path/ccdflist[0],ccdcal)

    oxyfull = [222,190]
    wh0 = (140,140)
    pstr=('0s','0.5s','1.0s','1.5s')
    lblring=array((89,87,85))
    ringmult=90-lblring

    for f,ring,ray,p in zip(ccdflist,rings,rays,pstr):
        imgfn = path/f
        outfn = path/('anno_' + f)
        fiducial(imgfn, a.xycrop[0], a.xycrop[1], outfn, ring, ray, p,
                 oxyfull,wh0,lblring,ringmult)

#%% sCMOS Figure 2
    flist = ('flame100615440.png','flame100615640.png','flame100615840.png',
             'flame100616040.png')
    #zero-based indexing, magnetic zenith
    oxyfull = [1328,876]
    pstr=('0s','0.2s','0.4s','0.6s') #times to annotate
    #half-width, half-height of the first oval
    # that is, 5 degreees elevation = wh0[1]
    wh0 = (252,234)
    lblring=array((88,86,84,82))
    axlim=(1000,1650,1200,550)
    ringmult=90-lblring

    if a.file: #user specified a single file
        imgfn = path/a.file
        outfn = path/('anno_' + a.file)
        fiducial(imgfn, a.xycrop[0], a.xycrop[1], outfn, rings[0], rays[0],
                 pstr[0],oxyfull,wh0,lblring,ringmult,axlim)
    else: #use default list
        for f,ring,ray,p in zip(flist,rings,rays,pstr):
            imgfn = path/f
            outfn = path/('anno_' + f)
            fiducial(imgfn, a.xycrop[0], a.xycrop[1], outfn, ring, ray,
                     p,oxyfull,wh0,lblring,ringmult,axlim)

    show()
