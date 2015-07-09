from __future__ import division,absolute_import
from numpy import logspace
import h5py
from os.path import join
from mpl_toolkits.mplot3d import Axes3D #needed for this file
#
from pymap3d.coordconv3d import ecef2aer, ecef2geodetic

def get1Dcut(cam,makeplot,progms,dbglvl):
    discardEdgepix = True #gets rid of duplicates beyond FOV of image that cause lsq estimation error
#%% determine slant range between other camera and magnetic zenith to evaluate at
    srpts = logspace(4.3,6.9,25) #4.5 had zero discards for hst0 #6.8 didn't quite get to zenith
#%% (0) load az/el data from Astrometry.net
    for c in cam:
        with h5py.File(cam[c].cal1Dfn,'r',libver='latest') as f:
            cam[c].doorient(f['/az'],f['/el'], f['/ra'], f['/dec'])
            #xind = fid['/x'].value;      #yind = fid['/y'].value
            #caltime = parser.parse(f['/timeFrame'])
        cam[c].toecef(srpts)

    #optional: plot ECEF of points between each camera and magnetic zenith (lying at az,el relative to each camera)
    plotLOSecef(cam,makeplot,progms,dbglvl)
#%% (2) get az,el of these points from camera to the other camera's points
    cam[0].az2pts,cam[0].el2pts,cam[0].r2pts = ecef2aer(cam[1].x2mz, cam[1].y2mz, cam[1].z2mz,
                                                             cam[0].lat, cam[0].lon, cam[0].alt_m)
    cam[1].az2pts,cam[1].el2pts,cam[1].r2pts = ecef2aer(cam[0].x2mz, cam[0].y2mz, cam[0].z2mz,
                                                             cam[1].lat, cam[1].lon, cam[1].alt_m)
#%% (3) find indices corresponding to these az,el in each image
        # and Least squares fit line to nearest points found in step 3
    for c in cam:
        cam[c].findClosestAzel(discardEdgepix)

#%%
    if dbglvl>2 and progms is not None:
        dbgfn = join(progms,'debugLSQ.h5')
        print('writing', dbgfn)
        with h5py.File(dbgfn,libver='latest') as fid:
            for c in cam:
                fid.create_dataset('/cam'+c+'/cutrow',data= cam[c].cutrow)
                fid.create_dataset('/cam'+c+'/cutcol', data = cam[c].cutcol)
                fid.create_dataset('/cam'+c+'/xpix', data = cam[c].xpix)
    return cam

def plotLOSecef(cam,makeplot,progms,dbglvl):
    from matplotlib.pyplot import figure
    if dbglvl>0:
        figecef = figure()
        clr = ['b','r','g','m']
        if dbglvl>1:
            import simplekml as skml
            kml1d = skml.Kml()


    for c in cam:
        if dbglvl>0: #SHOW PLOT
            axecef = figecef.gca(projection='3d')
            axecef.plot(xs=cam[c].x2mz, ys=cam[c].y2mz, zs=cam[c].z2mz, zdir='z',
                        color=clr[int(c)], label=c)
            axecef.set_title('LOS to magnetic zenith')

        if 'kmlell' in makeplot: #Write KML
            #convert LOS ECEF -> LLA
            loslat,loslon,losalt = ecef2geodetic(cam[c].x2mz,cam[c].y2mz,cam[c].z2mz)
            kclr = ['ff5c5ccd','ffff0000']
            #camera location points
            bpnt = kml1d.newpoint(name='HST'+c, description='camera ' +c + ' location',
                     coords=[(cam[c].lon,cam[c].lat)])
            bpnt.altitudemode = skml.AltitudeMode.clamptoground
            bpnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/paddle/pink-blank.png'
            bpnt.style.iconstyle.scale = 2.0
            #show cam to mag zenith los
            linestr = kml1d.newlinestring(name='')
            #TODO this is only first and last point without middle!
            linestr.coords = [(loslon[0],   loslat[0],  losalt[0]),
                              (loslon[-1], loslat[-1], losalt[-1])]
            linestr.altitudemode = skml.AltitudeMode.relativetoground
            linestr.style.linestyle.color = kclr[int(c)]

    if dbglvl>0:
        axecef.legend()
    if 'kmlell' in makeplot and progms is not None:
        kmlfn = join(progms,'debug1dcut.kmz')
        print('saving', kmlfn)
        kml1d.savekmz(kmlfn)


#
#def findClosestAzel(cam, discardEdgepix,dbglvl):
#    for c in cam:
#        azImg,elImg,azVec,elVec = cam[c].az, cam[c].el, cam[c].az2pts,cam[c].el2pts,
#
#        ny,nx = cam[c].ypix, cam[c].xpix
#
#        assert azImg.shape ==  elImg.shape
#        assert azVec.shape == elVec.shape
#        assert azImg.ndim == 2
#
#        npts = azVec.size  #numel
#        nearRow = empty(npts,dtype=int)
#        nearCol = empty(npts,dtype=int)
#        for ipt in range(npts):
#            #we do this point by point because we need to know the closest pixel for each point
#            errdist = absolute( hypot(azImg - azVec[ipt],
#                                       elImg - elVec[ipt]) )
#
## ********************************************
## THIS UNRAVEL_INDEX MUST BE ORDER = 'C'
#            nearRow[ipt],nearCol[ipt] = unravel_index(errdist.argmin(),(ny,nx),order='C')
##************************************************
#
#
#        if discardEdgepix:
#            edgeind = where(logical_or(logical_or(nearCol==0,nearCol == nx-1),
#                               logical_or(nearRow==0,nearRow==ny-1)) )[0]
#            nearRow = delete(nearRow,edgeind)
#            nearCol = delete(nearCol,edgeind)
#            if dbglvl>0: print('deleted',edgeind.size, 'edge pixels ')
#
#        cam[c].findLSQ(nearRow, nearCol)
#
#        if dbglvl>0:
#            clr = ['b','r','g','m']
#            ax = figure().gca()
#            ax.plot(nearCol,nearRow,color=clr[int(c)],label='cam'+c+'preLSQ',
#                    linestyle='None',marker='.')
#            ax.legend()
#            ax.set_xlabel('x'); ax.set_ylabel('y')
#            #ax.set_title('pixel indices (pre-least squares)')
#            ax.set_xlim([0,cam[c].az.shape[1]])
#            ax.set_ylim([0,cam[c].az.shape[0]])
#
#    return cam