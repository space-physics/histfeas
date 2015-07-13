#!/usr/bin/env python
# montage anno_flame10061*.png -trim -tile 4x1 -geometry +1+0  out.png
from matplotlib.patches import Ellipse
from matplotlib.pyplot import figure,show, Axes
import matplotlib
matplotlib.rcParams.update({'font.family':'sans-serif',
                           'font.sans-serif':'Arial',
                           'text.usetex':True})
from scipy.misc import imread
from numpy import array,cos,sin,radians
from os.path import expanduser, join

#half-width, half-height of the first oval
# that is, 5 degreees elevation = wh0[1]
wh0 = (252,234)

def main(imgfn,xcrop,ycrop,outfn,rings,rays,pstr):

    # geographic zenith, full frame video including axes/margins
    #%% setup rings
    xyr = []
    #zero-based indexing, geographic zenith
    #oxyfull = array([1242,297]) #by inspection of original 2013 Dahlgren video

    #zero-based indexing, magnetic zenith
    oxyfull = array([1328,876])
    #shift to account for post-processing crop of video (by me)
    oxy = array([oxyfull[0]-xcrop, oxyfull[1]-ycrop])

    xyr.append(array([wh0[0]*2/5-xcrop, wh0[1]*2/5-ycrop])) #by inspection

    #width,height of rest of ovals
    for i in range(2,5):
        xyr.append(i*xyr[0]+1) #multiple of radius
    #%% load image
    img = imread(imgfn)
    fg = figure()
    ax = Axes(fg, [0., 0., 1., 1.])
    ax. set_axis_off()
    fg.add_axes(ax)
    ax.axis((1000,1650,1200,550))
    ax.imshow(img)
    #ax.autoscale(True) #scale to image, not rings
    #%% plot rings
    if rings:
        ringincr = wh0[1]*2/5 #every 2 degrees
        for c in xyr:
            ax.add_patch(Ellipse(xy=oxy,
                             width=2*c[0], height=2*c[1],
                             fc='none',ec='white',lw=1))
        #%% annotate rings
        angtxt(ringincr, oxy,'88$^\circ$',ax)
        angtxt(ringincr*2,oxy,'86$^\circ$',ax)
        angtxt(ringincr*3,oxy,'84$^\circ$',ax)
        angtxt(ringincr*4,oxy,'82$^\circ$',ax)
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
    ax.text(1005,565,pstr,color='white',fontsize=50,
            va='top',ha='left',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))
    #%% finalize image
    #ax.axis('off')
    #fg.tight_layout

    #ax.set_ylim((500,1000))
    if outfn is not None:
        print('writing',outfn)
        fg.savefig(outfn,bbox_inches='tight',dpi=150)


def angtxt(radput,oxy,txt,ax):
    angput = radians(40)
    ax.text(wh0[0]/wh0[1]*radput*cos(angput)+oxy[0],
            radput*sin(angput)+oxy[1],
            txt,color='white',fontsize=45,
            va='top',ha='right',
            bbox=dict(boxstyle="round,pad=0.0",fc='black', alpha=0.25))

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='program to annotate cmos video avi')
    p.add_argument('path',help='path of images',type=str)
    p.add_argument('-c','--xycrop',help='origin of crop in x,y',nargs=2,default=(0,0),type=int)
    p.add_argument('-f','--file',help='name of single file to process',default=None,type=str)
    p.add_argument('-m','--magzen',help='azimuth from north (deg.) of magnetic zenith in image coordinate system',default=None,type=float)
    a = p.parse_args()

    flist = ('flame100615440.png','flame100615640.png','flame100615840.png',
             'flame100616040.png')

    pstr=('0s','0.2s','0.4s','0.6s') #times to annotate
    rings=(True,False,False,False)
    rays=(False,False,False,False) # set first True, and oxyfull to geographic zenith to plot ray pointing to magnetic zentih

    path = expanduser(a.path)

    if a.file is not None: #user specified a single file
        imgfn = join(path,a.file);
        outfn = join(path,'anno_' + a.file)
        main(imgfn, a.xycrop[0], a.xycrop[1], outfn, rings[0], rays[0], pstr[0])
    else: #use default list

        for f,ring,ray,p in zip(flist,rings,rays,pstr):
            imgfn = join(path,f);
            outfn = join(path,'anno_' + f)
            main(imgfn, a.xycrop[0], a.xycrop[1], outfn, ring, ray, p)

    show()
