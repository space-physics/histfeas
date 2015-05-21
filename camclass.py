from __future__ import print_function, division
from numpy import (linspace, fliplr, flipud, rot90, arange,
                   polyfit,polyval,rint,empty, isfinite,
                   absolute, hypot, logical_or, unravel_index, delete, where)
from os.path import expanduser
from dateutil.parser import parse
from scipy.signal import savgol_filter
from six import string_types
from numpy.random import poisson
from warnings import warn
#
from pymap3d.azel2radec import azel2radec
from pymap3d.haversine import angledist
from pymap3d.coordconv3d import aer2ecef


class Cam: #use this like an advanced version of Matlab struct
    def __init__(self,sim,cp,name,zmax,verbose):
        self.verbose = verbose
        self.name = name
#%%
        self.lat = cp['latWGS84']
        self.lon = cp['lonWGS84']
        self.alt_m = cp['Zkm']*1e3
        self.z_km = cp['Zkm']
        self.x_km = cp['Xkm']
        self.fn = cp['fn']
#        self.startTime = startTime
#        self.stopTime = stopTime
        self.nCutPix = int(cp['nCutPix'])
        self.xpix = cp['xPix']
        self.ypix = cp['yPix']
        self.xbin = cp['xBin']
        self.ybin = cp['yBin']
        self.Bincl = cp['Bincl']
        self.Bdecl = cp['Bdecl']
        self.Bepoch = cp['Bepoch'] #it's OK, I want to feedthru nan if xls cell is empty!

        if isinstance(self.Bepoch,string_types):
            self.Bepoch = parse(self.Bepoch)

        if isfinite(self.Bincl) and isfinite(self.Bdecl):
            self.Baz = 180. + self.Bdecl
            self.Bel = self.Bincl
        else:
            self.Baz = None; self.Bel = None

        self.rotCCW =    cp['rotCCW']
        self.transpose = cp['transpose'] == 1
        self.flipLR =    cp['flipLR'] == 1
        self.flipUD =   cp['flipUD'] == 1
        self.timeShiftSec = cp['timeShiftSec']

        self.plotminmax = (cp['plotMinVal'], cp['plotMaxVal'])
        self.fullstart = cp['fullFileStartUTC']

        self.intensityScaleFactor = cp['intensityScaleFactor']
        self.lowerthres = cp['lowerthres']

#%% check FOV and 1D cut sizes for sanity
        self.fovmaxlen = cp['FOVmaxLengthKM']

        if self.fovmaxlen > 10e3:
            print('sanityCheck: Your FOV length seems excessive > 10000 km')
        if self.nCutPix > 4096:
            print('sanityCheck: Program execution time may be excessive due to large number of camera pixels')
        if self.fovmaxlen < (1.5*zmax):
            print('sanityCheck: To avoid unexpected pixel/sky voxel intersection problems, make your candidate camera FOV at least 1.5 times longer than your maximum Z altitude.')

        self.boresightEl = cp['boresightElevDeg']
        self.arbfov = cp['FOVdeg']
#%%
        self.kineticSec = 1. / cp['frameRateHz']
#%% camera model
        """
        A model for sensor gain
        pedn is photoelectrons per data number
        This is used in fwd model and data inversion
        """
        self.pixarea_sqcm = cp['pixarea_sqcm']
        self.pedn = cp['pedn']
        self.ampgain = cp['ampgain']

        if isfinite(self.kineticSec) and isfinite(self.pixarea_sqcm) and isfinite(self.pedn):
            self.intens2dn = self.kineticSec * self.pixarea_sqcm * self.ampgain * self.ampgain / self.pedn
        else:
            self.intens2dn = 1
#%% sky mapping
        cal1Ddir = sim.cal1dpath
        cal1Dname = cp['cal1Dname']
        if isinstance(cal1Ddir,string_types) and isinstance(cal1Dname,string_types):
            self.cal1Dfn = expanduser(cal1Ddir + cal1Dname)

        self.raymap = sim.raymap
        if self.raymap == 'arbitrary':
            self.angle_deg = self.arbanglemap()
#        elif self.raymap == 'astrometry': #not OOPable at this time
#            self.angle_deg = self.astrometrymap()
#        else:
#            exit('*** unknown ray mapping method ' + self.raymap)
#%% pixel noise/bias
        self.noiselam = cp['noiseLam']
        self.ccdbias = cp['CCDBias']
        self.debiasData = cp['debiasData']
        self.smoothspan = cp['smoothspan']
        self.savgolOrder = cp['savgolOrder']


    def ingestcamparam(self,finf):
        self.SuperX=finf['superx']
        self.SuperY=finf['supery']
        self.Nmetadata = finf['nmetadata']
        self.BytesPerFrame= finf['bytesperframe']
        self.PixelsPerImage = finf['pixelsperimage']

    def arbanglemap(self):
        '''
        here's the center angle of each pixel ( + then - to go from biggest to
        smallest angle)  (e.g. 94.5 deg to 85.5 deg. -- sweeping left to right)
        '''
        #raySpacingDeg = self.arbfov / self.nCutPix
        maxAng = self.boresightEl + self.arbfov/2
        minAng = self.boresightEl - self.arbfov/2
        return linspace(maxAng, minAng, num=self.nCutPix, endpoint=True)

    def astrometrymap(self):
        pass
        # not sure if get1Dcut.py can be OOP'd
        #FOVangWidthDeg =pixAngleDeg[-1,iCam] - pixAngleDeg[0,iCam]
        #ModelRaySpacingDeg = np.mean( np.diff(pixAngleDeg[:,iCam],n=1) ) # for reference purposes

    def toecef(self,ranges):
        self.x2mz, self.y2mz, self.z2mz = aer2ecef(self.Baz,self.Bel,ranges,
                                                   self.lat,self.lon,self.alt_m)
    def doorient(self,az,el,ra,dec):
        if self.transpose:
            if self.verbose>0:
                print('tranposing cam #{} az/el/ra/dec data. '.format(self.name))
            az  = az.T
            el  = el.T
            ra  = ra.T
            dec = dec.T
        if self.flipLR:
            if self.verbose>0:
                print('flipping horizontally cam #{} az/el/ra/dec data.'.format(self.name))
            az  = fliplr(az)
            el  = fliplr(el)
            ra  = fliplr(ra)
            dec = fliplr(dec)
        if self.flipUD:
            if self.verbose>0:
                print('flipping vertically cam #{} az/el/ra/dec data.'.format(self.name))
            az  = flipud(az)
            el  = flipud(el)
            ra  = flipud(ra)
            dec = flipud(dec)
        if self.rotCCW != 0:
            if self.verbose>0:
                print('rotating cam #{} az/el/ra/dec data.'.format(self.name))
            az  = rot90(az, k = self.rotCCW)
            el  = rot90(el, k = self.rotCCW)
            ra  = rot90(ra, k = self.rotCCW)
            dec = rot90(dec,k = self.rotCCW)
        self.az = az
        self.el = el
        self.ra = ra
        self.dec = dec

    def debias(self,data):
        if isfinite(self.debiasData):
            if self.verbose>0:
                print('Debiasing Data for Camera #{}'.format(self.name) )
            data -= self.debiasData
        return data

    def donoise(self,data):
         noisy = data.copy()

         if isfinite(self.noiselam):
             if self.verbose>0:
                 print('adding Poisson noise with \lambda={:0.1f} to camera #{}'.format(self.noiselam,self.name))
             dnoise = poisson(lam=self.noiselam,size=self.nCutPix)
             noisy += dnoise
         else:
             dnoise = None


         if isfinite(self.ccdbias):
             if self.verbose>0:
                 print('adding bias {:0.1e} to camera #{}'.format(self.ccdbias,self.name))
             noisy += self.ccdbias

         # kept for diagnostic purposes
         self.raw = data #these are untouched pixel intensities
         self.dnoise = dnoise
         self.noisy = noisy

         return noisy

    def dosmooth(self,data):
        if self.smoothspan > 0 and self.savgolOrder>0:
            if self.verbose>0:
                print('Smoothing Data for Camera #{}'.format(self.name))
            data= savgol_filter(data, self.smoothspan, self.savgolOrder)
        return data #LEAVE THIS AS SEPARATE LINE!

    def fixnegval(self,data):
        negDataInd = data<0
        if (self.verbose and negDataInd.any()) or negDataInd.sum()>0.1*self.nCutPix:
            warn('Setting {} negative Data values to 0 for Camera #{}'.format(negDataInd.sum(), self.name))
        
        data[negDataInd] = 0

        self.nonneg = data
        return data
    def scaleintens(self,data):
        if isfinite(self.intensityScaleFactor) and self.intensityScaleFactor !=1 :
            if self.verbose>0:
                print('Scaling data to Cam #' + self.name + ' by factor of ',self.intensityScaleFactor )
            data *= self.intensityScaleFactor
            #assert isnan(data).any() == False
        return data

    def dolowerthres(self,data):
        if isfinite(self.lowerthres):
            print('Thresholding Data to 0 when DN <',self.lowerthres, 'for Camera #' + self.name )
            data[ data < self.lowerthres ] = 0
        return data

    def findLSQ(self,nearrow,nearcol):
        polycoeff = polyfit(nearcol,nearrow,deg=1,full=False)
        #columns (x)  to cut from picture
        cutcol = arange(self.xpix,dtype=int) #not range
        #rows (y) to cut from picture
        cutrow = rint(polyval(polycoeff,cutcol)).astype(int)
        assert (cutrow>=0).all() and (cutrow<self.ypix).all()

        #angle from magnetic zenith corresponding to those pixels
        rapix =  self.ra[cutrow, cutcol]
        decpix = self.dec[cutrow, cutcol]
        raMagzen,decMagzen = azel2radec(self.Baz,self.Bel,self.lat,self.lon,self.Bepoch)
        if self.verbose>0: print('mag. zen. ra/dec {} {}'.format(raMagzen,decMagzen))

        angledist_deg = angledist(raMagzen,decMagzen,rapix,decpix)
        # put distances into a 90-degree fan beam
        angle_deg = empty(self.xpix,dtype=float)
        MagZenInd = angledist_deg.argmin() # whether slightly positive or slightly negative, this should be OK

        angle_deg[MagZenInd:] = 90. + angledist_deg[MagZenInd:]
        angle_deg[:MagZenInd] = 90. - angledist_deg[:MagZenInd]

        self.angle_deg = angle_deg
        self.angleMagzenind = MagZenInd
        self.cutrow = cutrow
        self.cutcol = cutcol

        if self.verbose>0:
            from matplotlib.pyplot import figure
            clr = ['b','r','g','m']
            ax=figure().gca()
            ax.plot(cutcol,cutrow,color=clr[self.name],
                    label=self.name, linestyle='-')
            ax.legend()
            ax.set_xlabel('x'); ax.set_ylabel('y')
            ax.set_title('polyfit with computed ray points')

            ax =figure().gca()
            ax.plot(angle_deg,color=clr[self.name],label='cam'+self.name,
                    linestyle='None',marker='.')
            ax.legend()
            ax.set_xlabel('x-pixel'); ax.set_ylabel('$\theta$ [deg.]')
            ax.set_title('angle from magnetic zenith $\theta$')

    def findClosestAzel(self,discardEdgepix):
        assert self.az.shape ==  self.el.shape
        assert self.az2pts.shape == self.el2pts.shape
        assert self.az.ndim == 2

        npts = self.az2pts.size  #numel
        nearRow = empty(npts,dtype=int)
        nearCol = empty(npts,dtype=int)
        #FIXME consider scipy.spatial.distance.cdist()
        for ipt in range(npts):
            #we do this point by point because we need to know the closest pixel for each point
            errdist = absolute( hypot(self.az - self.az2pts[ipt],
                                      self.el - self.el2pts[ipt]) )

    # ********************************************
    # THIS UNRAVEL_INDEX MUST BE ORDER = 'C'
            nearRow[ipt],nearCol[ipt] = unravel_index(errdist.argmin(),
                                                      (self.ypix, self.xpix),order='C')
    #************************************************


        if discardEdgepix:
            edgeind = where(logical_or(logical_or(nearCol==0,nearCol == self.xpix-1),
                                       logical_or(nearRow==0,nearRow == self.ypix-1)) )[0]
            nearRow = delete(nearRow,edgeind)
            nearCol = delete(nearCol,edgeind)
            if self.verbose>0: print('deleted',edgeind.size, 'edge pixels ')

        self.findLSQ(nearRow, nearCol)

        if self.verbose>0:
            from matplotlib.pyplot import figure
            clr = ['b','r','g','m']
            ax = figure().gca()
            ax.plot(nearCol,nearRow,color=clr[int(self.name)],label='cam{}preLSQ'.format(self.name),
                    linestyle='None',marker='.')
            ax.legend()
            ax.set_xlabel('x'); ax.set_ylabel('y')
            #ax.set_title('pixel indices (pre-least squares)')
            ax.set_xlim([0,self.az.shape[1]])
            ax.set_ylim([0,self.az.shape[0]])