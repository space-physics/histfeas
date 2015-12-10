from __future__ import print_function, division,absolute_import
import h5py
#from numba import jit
#from numbapro import vectorize
from numpy import empty,ones,ravel_multi_index,hypot,zeros,in1d,array
from scipy.sparse import dok_matrix,issparse
from shutil import copy2
# local
from cvutils.lineClipping import cohensutherland

'''
 Michael Hirsch

 This code is a nest of four for loops
 We want to implement b = Lv, where
   b is an nCam*nPixel column-major vector of observed brightness (data numbers)
   L is an nCam*nPixel x sx*sz dimension projection matrix
   v is an sx*sz column-major vector of auroral VER
 by building up b,L,v in this way, the tomographic forward projection is implemented
 in Python numpy by b = L.dot(v), which in Matlab is b = L*v

 let's start by describing b,L,v.
 We have made a mesh in the sky, let's say 41 elements from -10 to 10km with 0.5km horizontal spacing (x)
 and say 16 elements from 90 to 400 km with 20km vertical spacing (z)
 That means v will be a 16*41 = 656 element vector and L will have 656 columns
 now let's say we have 2 cameras (nCam=2) each with 512 pixels in a 1-D cut of the sky through a common volume (so carefully selected).
 then L will have 2*512 = 1024 rows, thus L is a 1024x656 matrix.
 finally, b will be 2*512 = 1024 element vector.

 next, the for loops from the inside out
 The innermost loop (z) goes over each altitude from bottom to top (90km to 390km in our example)
 the next innermost loop (x) goes over each horizontal bin from left to right (-10km to 10km in our example)
 the next loop (k) is for each pixel of the particular camera
 the outer loop (cam) is for each camera

 in words, this program implements the Cohen-Sutherland line-clipping algorithm over each of the sky
 pixels from bottom to top, left to right, for each pixel of each camera. it finds the
 line length "ell" for each element of L. This is how we implement the line integral.
'''

def EllLineLength(Fwd,xFOVpixelEnds,zFOVpixelEnds,Np,sim,makePlots,dbglvl):
    plotEachRay=False
    writeRays = False #write pixel rays to hdf5 for viewing

#%% convenience variables (easier to type)
    sx =  Fwd['sx']
    sz =  Fwd['sz']
    nCam = sim.nCamUsed
    xpc = Fwd['xPixCorn']
    zpc = Fwd['zPixCorn']
    maxNell = Fwd['maxNell']
    assert xpc.size == sx + 1
    assert zpc.size == sz + 1
#%% preallocation
#This goes OUTSIDE all loops!
    if dbglvl>0:
        print('EllLineLength: Number of 1D pixels in cut: ',Np)
        print('x: ', end=''); print(Fwd['x'])
        print('xpc: ',end=''); print(xpc)
        print('z: ', end=''); print(Fwd['z'])
        print('zpc: ',end=''); print(zpc)

    xzplot =None
    L = goCalcEll(maxNell,nCam,Np,sz,sx,xpc,zpc,xFOVpixelEnds,zFOVpixelEnds,
                          sim.allCamXkm,sim.allCamZkm,plotEachRay,dbglvl)

    #%% write results to HDF5 file
    if sim.savefwdL:
        doSaveEll(L,Fwd,sim,xFOVpixelEnds,zFOVpixelEnds,writeRays)

    if 'ell' in makePlots and plotEachRay and xzplot:
        plotEll(nCam,xFOVpixelEnds,zFOVpixelEnds,sim.xCam,sim.zCam,Np,xpc,zpc,
                sz,sx,xzplot,sim.FwdLfn,plotEachRay,makePlots,
                (None,None,None,None,None,None))
    if issparse(L):
        L = L.tocsc()
    return L

def goCalcEll(maxNell,nCam,Np,sz,sx,xpc,zpc,xFOVpixelEnds,zFOVpixelEnds,
                xCam,zCam,plotEachRay,dbglvl=0):
    usesparse=True
    if usesparse:
        L = dok_matrix(( Np*nCam,sz*sx),dtype=float) #sparse
    else:
        L = zeros( ( Np*nCam,sz*sx),dtype=float ,order='F') #dense


    print('Dimensions of L:',L.shape,' sz=',sz, '  sx=',sx )

    Lcol =   empty(maxNell, dtype=int) #we'll truncate this to actual length at end
    tmpEll = empty(maxNell, dtype=float) #we'll truncate this to actual length at end
    xzplot = [] #we'll append to this

    L = loopEll(Np,sz,sx,xpc,zpc,xFOVpixelEnds,zFOVpixelEnds,
                             xCam,zCam,nCam,plotEachRay,
                             L,Lcol,tmpEll,xzplot) #numba

    if usesparse:
        return L.tocsc()
    else:
        return L

#@jit(['float64[:,:](int64,int64,int64,float64[:],float64[:],float64[:,:],float64[:,:],float64[:],float64[:],int64[:],bool_,float64[:,:],int64[:],float64[:],float64[:])'])
def loopEll(Np,sz,sx,xpc,zpc,xFOVpixelEnds,zFOVpixelEnds,
             xCam,zCam,nCam,plotEachRay,L,Lcol,tmpEll,xzplot):
#%%let's compute intersections!!
    #FIXME assumes all cameras have same # of pixels
    inttot = 0

    ''' We MUST have all cameras enabled for this computation (as verified in observeVolume)'''
#    try:
    for iCam in range(nCam):
        xfov = xFOVpixelEnds[:,iCam]; zfov = zFOVpixelEnds[:,iCam]
        for k in range(Np): #FIXME assumes all cam same 1D cut pixel length
            nHitsThisPixelRay = 0
            for xInd in range(sx):
                for zInd in range(sz):
                    '''get the vertices of the intersection between voxel and pixel
                    for the kth pixel FOV, find the intersection with each sky polygon
                     '''
                    x1, y1, x2, y2 = cohensutherland(
                             xpc[xInd], zpc[zInd + 1],
                              xpc[xInd + 1], zpc[zInd],
                              xCam[iCam], zCam[iCam],
                              xfov[k], zfov[k])
                    if x1 is not None:
                        #assert nHitsThisPixelRay <= maxNell #removed for performance
                        Lcol[nHitsThisPixelRay] = ravel_multi_index((zInd,xInd),dims=(sz,sx),order='F')
                        tmpEll[nHitsThisPixelRay] = hypot(x2-x1,y2-y1)
                        nHitsThisPixelRay += 1
                        if plotEachRay:
                            xzplot.append([x1,x2,y1,y2])
            if nHitsThisPixelRay > 0: #this is under "for k" level
                inttot += nHitsThisPixelRay
                L[iCam * Np + k, Lcol[:nHitsThisPixelRay]] = tmpEll[:nHitsThisPixelRay]
            if k%100 == 0: #arbitrary display update interval #this is under "for k" level
                print('Camera #{}: {:0.0f}% complete, found {} intersections.'.format(iCam,1.0*k/Np*100,inttot))

#    except IndexError:
#        exit('**** ERROR on iCam='+str(iCam)+', k=' +str(k)+', xInd='+str(xInd)+', zInd='+str(zInd)+' *****')
    print('Total number of intersections found:',inttot)
    return L#,xzplot

def doSaveEll(L,Fwd,sim,xFOVpixelEnds,zFOVpixelEnds,writeRays):
    print('writing {}'.format(sim.FwdLfn))
    if issparse(L):
        L = L.todense()
    with h5py.File(str(sim.FwdLfn),'w',libver='latest') as fid:
        h5L = fid.create_dataset("/L",data=L,compression="gzip");  h5L.attrs['Units'] = 'kilometers'
        h5Fwdx = fid.create_dataset("/Fwd/x",data=Fwd['x']); h5Fwdx.attrs['Units'] = 'kilometers'
        h5Fwdz = fid.create_dataset("/Fwd/z",data=Fwd['z']); h5Fwdz.attrs['Units'] = 'kilometers'
        h5FwdxPC = fid.create_dataset("/Fwd/xPixCorn",data=Fwd['xPixCorn']); h5FwdxPC.attrs['Units'] = 'kilometers'
        h5FwdzPC = fid.create_dataset("/Fwd/zPixCorn",data=Fwd['zPixCorn']); h5FwdzPC.attrs['Units'] = 'kilometers'

        #h5ObsPA =  fid.create_dataset("/Obs/pixAngle",data=pixAngleDeg);  h5ObsPA.attrs['Units'] = 'Degrees'
        h5ObsxFPE = fid.create_dataset("/Obs/xFOVpixelEnds",data=xFOVpixelEnds); h5ObsxFPE.attrs['Units'] = 'kilometers'
        h5ObszFPE = fid.create_dataset("/Obs/zFOVpixelEnds",data=zFOVpixelEnds); h5ObszFPE.attrs['Units'] = 'kilometers'
        h5xCam = fid.create_dataset('/Obs/xCam',data=sim.allCamXkm); h5xCam.attrs['Units'] = 'kilometers'
        h5zCam = fid.create_dataset('/Obs/zCam',data=sim.allCamZkm); h5zCam.attrs['Units'] = 'kilometers'
    copy2(str(sim.FwdLfn), sim.cal1dpath)


def plotEll(nCam,xFOVpixelEnds,zFOVpixelEnds,xCam,zCam,Np,xpc,zpc,sz,sx,
            xzplot,EllFN,plotEachRay,makeplot,vlim):

    from matplotlib.pyplot import figure, draw, pause
    from matplotlib.ticker import MultipleLocator
    decimfactor = 8 #plot every Nth ray
    clrs = ['r','g','y','m']
    afs = None#20
    tkfs = None#20
    tfs = None#22

    fg = figure()
    ax = fg.gca()
    #this pcolormesh shows the model grid, with a white background per JLS
    ax.pcolormesh(xpc,zpc,ones((sz,sx)), cmap='bone', vmin= 0, vmax=1,
                  edgecolor='k', linewidth=0.001)


    if plotEachRay:
        print('plotting viewing geometry')
        for ixz in range(len(xzplot)):
            ax.plot(x=xzplot[ixz][:2],
                    y=xzplot[ixz][-2:],
                    color='red')
            draw() #need plt.pause(0.01) just after this
            pause(0.01) #need plt.draw() just before this
    else:
        for iCam in range(nCam):
            for iray in range(0,Np,decimfactor):
                #DO NOT SPECIFY x=... y=... or NO LINES appear!!
                ax.plot([xCam[iCam],xFOVpixelEnds[iray,iCam]],
                        [zCam[iCam],zFOVpixelEnds[iray,iCam]],
                        color=clrs[iCam])
    if vlim[0] is None:
        ax.set_xlim((xpc[0], xpc[-1]))
        ax.set_ylim((0,      zpc[-1]))
    else:
        ax.set_xlim(vlim[:2])
        ax.set_ylim(vlim[2:4])
    ax.set_xlabel('$B_\perp$ [km]', fontsize = afs)
    ax.set_ylabel('$B_\parallel$ [km]', fontsize = afs)
    ax.set_title('Viewing geometry for $L$', fontsize = tfs)# + EllFN)

    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax.tick_params(axis='both', which='both', labelsize=tkfs, direction='out')

#%% write plot
    tmpl = array(('eps','jpg','png','pdf'))
    used = in1d(tmpl,makeplot)
    if used.any():
        figpng = EllFN + '.' + tmpl[used][0]
        print('saving', figpng)
        fg.savefig(figpng,bbox_inches='tight',dpi=600)  # this is slow and async
