from pathlib import Path
import logging
from numpy import (linspace, fliplr, flipud, rot90, arange,
                   polyfit,polyval,rint,empty, isfinite, isclose,
                   absolute, hypot, unravel_index,logical_not,array)
from datetime import datetime
from pytz import UTC
from dateutil.parser import parse
from scipy.signal import savgol_filter
from numpy.random import poisson
import h5py
from astropy.coordinates.angle_utilities import angular_separation
from astropy import units as u
#
from pymap3d.azel2radec import azel2radec
from pymap3d.coordconv3d import aer2ecef
from dascutils.readDASCfits import readDASC
from themisasi.fov import mergefov

epoch = datetime(1970,1,1,tzinfo=UTC)

verbose=0


class Cam: #use this like an advanced version of Matlab struct
    def __init__(self,sim,cp,name,zmax,makeplot,verbose=0):
        self.verbose = verbose

        self.usecam = bool(cp['useCam'])
        if not self.usecam and sim.realdata and name.lower() == 'asi':
            self.fn = list(Path(cp['fn']).expanduser().glob('*.FITS'))
            self.name='asi'

            self.clim = [None]*2
            if isfinite(cp['plotMinVal']): self.clim[0] =  cp['plotMinVal']
            if isfinite(cp['plotMaxVal']): self.clim[1] =  cp['plotMaxVal']

            _,(self.az,self.el),self.lla,_ = readDASC(self.fn[0],cp['azcalfn'],cp['elcalfn'])


            if 'realvid' in makeplot:
                if 'h5' in makeplot:
                    print('fov: writing {}'.format(sim.fovfn))
                    self.hlrows,self.hlcols = mergefov(None,self.lla,self.az,self.el,None,None,
                                   ['../histutils/cal/hst0cal.h5','../histutils/cal/hst1cal.h5'],
                                   projalt=110e3,site='DASC')
                    with h5py.File(sim.fovfn,'w',libver='latest') as H:
                        H['/rows'] = array(self.hlrows) # Ncam x 4 x Nx  (list,list,ndarray)
                        H['/cols'] = array(self.hlcols)
                else:
                    with h5py.File(sim.fovfn,'r',libver='latest') as H:
                        self.hlrows = H['/rows'].value
                        self.hlcols = H['/cols'].value

            return

        elif self.usecam:
            self.name = int(name)
#%%
        self.nCutPix = int(cp['nCutPix'])

        self.Bincl = cp['Bincl']
        self.Bdecl = cp['Bdecl']
        self.Bepoch = cp['Bepoch'] #it's OK, I want to feedthru nan if xls cell is empty!

        if isinstance(self.Bepoch,str):
            self.Bepoch = parse(self.Bepoch)

        if isfinite(self.Bincl) and isfinite(self.Bdecl):
            self.Baz = 180. + self.Bdecl
            self.Bel = self.Bincl
        else:
            self.Baz = None; self.Bel = None

        self.timeShiftSec = cp['timeShiftSec'] if isfinite(cp['timeShiftSec']) else 0.

        self.clim = [None]*2
        if isfinite(cp['plotMinVal']): self.clim[0] =  cp['plotMinVal']
        if isfinite(cp['plotMaxVal']): self.clim[1] =  cp['plotMaxVal']

        self.intensityScaleFactor = cp['intensityScaleFactor']
        self.lowerthres = cp['lowerthres']

#%% check FOV and 1D cut sizes for sanity
        self.fovmaxlen = cp['FOVmaxLengthKM']

        if self.fovmaxlen > 10e3:
            logging.warning('sanityCheck: Your FOV length seems excessive > 10000 km')
        if self.nCutPix > 4096:
            logging.warning('sanityCheck: Program execution time may be excessive due to large number of camera pixels')
        if self.fovmaxlen < (1.5*zmax):
            logging.warning('sanityCheck: To avoid unexpected pixel/sky voxel intersection problems, make your candidate camera FOV at least 1.5 times longer than your maximum Z altitude.')

        self.boresightEl = cp['boresightElevDeg']
        self.arbfov = cp['FOVdeg']

#%% sky mapping
        cal1Ddir = sim.rootdir/sim.cal1dpath
        cal1Dname = cp['cal1Dname']
        if isinstance(cal1Ddir,Path) and isinstance(cal1Dname,str):
            self.cal1Dfn = (cal1Ddir / cal1Dname).expanduser()

        self.raymap = sim.raymap
        if not sim.realdata and self.raymap == 'arbitrary': #don't use realdata with arbitrary, doesn't make sense and causes flipped brightness data when loading--you've been warned!
            self.angle_deg = self.arbanglemap()
#        elif self.raymap == 'astrometry': #not OOPable at this time
#            self.angle_deg = self.astrometrymap()
#        else:
#            exit('*** unknown ray mapping method ' + self.raymap)
#%% pixel noise/bias
        try:
            self.noiselam = cp['noiseLam']
            self.ccdbias = cp['CCDBias']
        except KeyError: #realdata
            pass

        self.debiasData = cp['debiasData']
        self.smoothspan = cp['smoothspan']
        self.savgolOrder = cp['savgolOrder']

        """ expects an HDF5 .h5 file"""

        # data file name
        if sim.realdata:
            self.fn = (sim.realdatapath / cp['fn']).expanduser()

            with h5py.File(str(self.fn),'r',libver='latest') as f:
                self.filestartutc = f['/ut1_unix'][0]
                self.filestoputc  = f['/ut1_unix'][-1]
                self.ut1unix      = f['/ut1_unix'].value + self.timeShiftSec
                self.supery,self.superx = f['/rawimg'].shape[1:]

                p = f['/params']
                self.kineticsec   = p['kineticsec']
                self.rotccw       = p['rotccw']
                self.transpose    = p['transpose'] == 1
                self.fliplr       = p['fliplr'] == 1
                self.flipud       = p['flipud'] == 1

                c = f['/sensorloc']
                self.lat   = c['lat'][0]
                self.lon   = c['lon'][0]
                self.alt_m = c['alt_m'][0]
        else: #sim
            self.kineticsec = cp['kineticsec'] #simulation
            self.alt_m = cp['Zkm']*1000
            self.x_km = cp['Xkm']

#%% camera model
        """
        A model for sensor gain
        pedn is photoelectrons per data number
        This is used in fwd model and data inversion
        """
        hasgainparam = isfinite(self.kineticsec) and isfinite(cp['pixarea_sqcm']) and isfinite(cp['pedn']) and isfinite(cp['ampgain'])

        if hasgainparam:
            self.dn2intens = cp['pedn'] / (self.kineticsec * cp['pixarea_sqcm'] * cp['ampgain'])
            if sim.realdata:
                self.intens2dn = 1
            else:
                self.intens2dn = 1/self.dn2intens
        else: #this will give tiny ver and flux
            self.intens2dn = self.dn2intens = 1

#%% summary
        if sim.realdata:
            logging.info('cam{} timeshift: {} seconds'.format(self.name,self.timeShiftSec))

            logging.info('Camera {} start/stop UTC: {} / {}, {} frames.'.format(
                                                              self.name,
                                                              self.filestartutc,
                                                              self.filestoputc,
                                                              self.ut1unix.size))


    def arbanglemap(self):
        '''
        here's the center angle of each pixel ( + then - to go from biggest to
        smallest angle)  (e.g. 94.5 deg to 85.5 deg. -- sweeping left to right)
        '''
        #raySpacingDeg = self.arbfov / self.nCutPix
        maxAng = self.boresightEl + self.arbfov/2
        minAng = self.boresightEl - self.arbfov/2
        angles=linspace(maxAng, minAng, num=self.nCutPix, endpoint=True)
        assert isclose(angles[0],maxAng) & isclose(angles[-1],minAng)
        return angles

    def astrometrymap(self):
        pass
        # not sure if get1Dcut.py can be OOP'd
        #FOVangWidthDeg =pixAngleDeg[-1,iCam] - pixAngleDeg[0,iCam]
        #ModelRaySpacingDeg = np.mean( np.diff(pixAngleDeg[:,iCam],n=1) ) # for reference purposes

    def toecef(self,ranges):
        self.x2mz, self.y2mz, self.z2mz = aer2ecef(self.Baz,self.Bel,ranges,
                                                   self.lat,self.lon,self.alt_m)

    #TODO put doorientimage and doorient in one function for safety
    def doorientimage(self,frame):
        if self.transpose:
            if frame.ndim==3:
                frame = frame.transpose(0,2,1)
            elif frame.ndim==2:
                frame = frame.T
            else:
                raise ValueError('ndim==2 or 3')
        # rotate -- note if you use origin='lower', rotCCW -> rotCW !
         #rotate works with first two axes
        if self.rotccw: #NOT isinstance integer_types!
            if frame.ndim == 3:
                frame = rot90(frame.transpose(1,2,0),k=self.rotccw).transpose(2,0,1)
            elif frame.ndim == 2:
                frame = rot90(frame,k=self.rotccw)
            else:
                raise ValueError('ndim==2 or 3')
        # flip
        if self.fliplr:
            if frame.ndim == 3:
                frame = fliplr(frame.transpose(1,2,0)).transpose(2,0,1)
            elif frame.ndim==2:
                frame = fliplr(frame)
            else:
                raise ValueError('ndim==2 or 3')
        if self.flipud:
            if frame.ndim == 3:
                frame = flipud(frame.transpose(1,2,0)).transpose(2,0,1)
            elif frame.ndim==2:
                frame = flipud(frame)
            else:
                raise ValueError('ndim==2 or 3')
        return frame

    def doorient(self,az,el,ra,dec):
        if self.transpose:
            logging.debug('tranposing cam #{} az/el/ra/dec data. '.format(self.name))
            az  = az.T
            el  = el.T
            ra  = ra.T
            dec = dec.T
        if self.fliplr:
            logging.debug('flipping horizontally cam #{} az/el/ra/dec data.'.format(self.name))
            az  = fliplr(az)
            el  = fliplr(el)
            ra  = fliplr(ra)
            dec = fliplr(dec)
        if self.flipud:
            logging.debug('flipping vertically cam #{} az/el/ra/dec data.'.format(self.name))
            az  = flipud(az)
            el  = flipud(el)
            ra  = flipud(ra)
            dec = flipud(dec)
        if self.rotccw != 0:
            logging.debug('rotating cam #{} az/el/ra/dec data.'.format(self.name))
            az  = rot90(az, k = self.rotccw)
            el  = rot90(el, k = self.rotccw)
            ra  = rot90(ra, k = self.rotccw)
            dec = rot90(dec,k = self.rotccw)
        self.az = az
        self.el = el
        self.ra = ra
        self.dec = dec

    def debias(self,data):
        if hasattr(self,'debiasData') and isfinite(self.debiasData):
            logging.debug('Debiasing Data for Camera #{} by -{}'.format(self.name,self.debiasData) )
            data -= self.debiasData
        return data

    def donoise(self,data):
         noisy = data.copy()

         if hasattr(self,'noiselam') and isfinite(self.noiselam):
             logging.info('adding Poisson noise with \lambda={:0.1f} to camera #{}'.format(self.noiselam,self.name))
             dnoise = poisson(lam=self.noiselam,size=self.nCutPix)
             noisy += dnoise
         else:
             dnoise = None


         if hasattr(self,'ccdbias') and isfinite(self.ccdbias):
             logging.info('adding bias {:0.1e} to camera #{}'.format(self.ccdbias,self.name))
             noisy += self.ccdbias

         # kept for diagnostic purposes
#         self.raw = data #these are untouched pixel intensities
         self.dnoise = dnoise
         self.noisy = noisy

         return noisy

    def dosmooth(self,data):
        assert isfinite(data).all(),'NaN leaked into brightness data, savgol cannot handle NaN'
        if self.smoothspan > 0 and self.savgolOrder>0:
            logging.debug('Smoothing Data for Camera #{}'.format(self.name))
            data= savgol_filter(data, self.smoothspan, self.savgolOrder)
        return data #LEAVE THIS AS SEPARATE LINE!

    def fixnegval(self,data):
        mask = data<0
        if mask.sum()>0.2*self.nCutPix:
            logging.info('Setting {} negative Data values to 0 for Camera #{}'.format(mask.sum(), self.name))

        data[mask] = 0

#        self.nonneg = data
        return data

    def scaleintens(self,data):
        if isfinite(self.intensityScaleFactor) and self.intensityScaleFactor !=1 :
            logging.info('Scaling data to Cam #{} by factor of {}'.format(self.name,self.intensityScaleFactor))
            data *= self.intensityScaleFactor
            #assert isnan(data).any() == False
        return data

    def dolowerthres(self,data):
        if isfinite(self.lowerthres):
            print('Thresholding Data to 0 when DN < {} for Camera {}'.format(self.lowerthres,self.name))
            data[ data < self.lowerthres ] = 0
        return data

    def findLSQ(self,nearrow,nearcol):
        polycoeff = polyfit(nearcol,nearrow,deg=1,full=False)
        #columns (x)  to cut from picture
        cutcol = arange(self.superx,dtype=int) #not range
        #rows (y) to cut from picture
        cutrow = rint(polyval(polycoeff,cutcol)).astype(int)
        assert (cutrow>=0).all() and (cutrow<self.supery).all(),'impossible least squares fit for 1-D cut\n is your video orientation correct? check the params of video hdf5 file'
        # DONT DO THIS: cutrow.clip(0,self.supery,cutrow)

        #angle from magnetic zenith corresponding to those pixels
        rapix =  self.ra[cutrow, cutcol]
        decpix = self.dec[cutrow, cutcol]
        raMagzen,decMagzen = azel2radec(self.Baz,self.Bel,self.lat,self.lon,self.Bepoch)
        logging.info('mag. zen. ra/dec {} {}'.format(raMagzen,decMagzen))

        angledist = angular_separation(raMagzen*u.deg,decMagzen*u.deg,rapix*u.deg,decpix*u.deg)
        angledist = angledist.to(u.deg).value
        # put distances into a 90-degree fan beam
        angle_deg = empty(self.superx,dtype=float)
        MagZenInd = angledist.argmin() # whether minimum angle distance from MZ is slightly positive or slightly negative, this should be OK

        angle_deg[MagZenInd:] = 90. + angledist[MagZenInd:]
        angle_deg[:MagZenInd] = 90. - angledist[:MagZenInd]

        self.angle_deg = angle_deg
        self.angleMagzenind = MagZenInd
        self.cutrow = cutrow
        self.cutcol = cutcol

        if verbose>0:
            from matplotlib.pyplot import figure
            clr = ['b','r','g','m']
            ax=figure().gca()
            ax.plot(cutcol,cutrow,color=clr[self.name],
                    label=self.name, linestyle='-')
            ax.legend()
            ax.set_xlabel('x'); ax.set_ylabel('y')
            ax.set_title('polyfit with computed ray points')

            ax =figure().gca()
            ax.plot(angle_deg,color=clr[self.name],label='cam {}'.format(self.name),
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
        # can be FAR FAR faster than scipy.spatial.distance.cdist()
        for i in range(npts):
            #we do this point by point because we need to know the closest pixel for each point
            errdist = absolute( hypot(self.az - self.az2pts[i],
                                      self.el - self.el2pts[i]) )

    # ********************************************
    # THIS UNRAVEL_INDEX MUST BE ORDER = 'C'
            nearRow[i],nearCol[i] = unravel_index(errdist.argmin(), self.az.shape, order='C')
    #************************************************


        if discardEdgepix:
            mask = logical_not(((nearCol==0) | (nearCol == self.az.shape[1]-1)) |
                               ((nearRow==0) | (nearRow == self.az.shape[0]-1)))
            nearRow = nearRow[mask]
            nearCol = nearCol[mask]

            self.findLSQ(nearRow, nearCol)

        if verbose>0:
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
