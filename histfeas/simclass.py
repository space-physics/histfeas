from pathlib import Path
import logging
from hashlib import md5
from numpy import asarray,where,arange,isfinite,ceil,hypot,atleast_1d,fromstring
import numpy as np # needed for all
from datetime import datetime
from dateutil.parser import parse
#
from transcarread.readionoinit import getaltgrid
from .camclass import splitconf

class Sim:

    def __init__(self,sp,ap,ntimeslice,overrides,makeplot,odir):
#%% root directory not here
        try:
            self.rootdir = Path(overrides['rootdir']).expanduser()
        except KeyError:
            self.rootdir = Path('')
#%% how many cameras in use, and which ones?
        self.camnames = sp['cam']['name'].split(',')
        self.useCamBool = fromstring(sp['cam']['useCam'],dtype=bool,sep=',')

        try:
            usecamreq = asarray(overrides['cam'])
            if usecamreq[0] is not None: #override spreadsheet
                logging.warning(' Overriding INI, using cameras: {}'.format(usecamreq))
                for i,civ in enumerate(self.useCamBool): # this might be a silly indexing method, but works
                    if (i==usecamreq).any():
                        self.useCamBool[i] = True
                    else:
                        self.useCamBool[i] = False
            if usecamreq[0] is not None:
                assert np.all(where(self.useCamBool)[0] == usecamreq) #not .all() in case of different length
        except KeyError:
            pass #normal case
#%% camera position override check (work done in sanityCheck.setupCam)
        try:
            self.camxreq = overrides['camx']
            if self.camxreq[0] and len(self.camxreq) != self.nCamUsed:
                raise ValueError('must specify same number of x-loc and used cameras')
        except (TypeError,KeyError): #cam override not specified
            self.camxreq = [None]
#%%
        self.nCamUsed = self.useCamBool.sum() #result is an int

        nCutPix = fromstring(sp['cam']['nCutPix'],dtype=int,sep=',') #FIXME someday allow different # of pixels..
        assert (nCutPix[self.useCamBool] == nCutPix[self.useCamBool][0]).all(),'all cameras must have same 1D cut length'
        self.ncutpix = nCutPix[self.useCamBool][0]

        if 'realvid' in makeplot:
            self.fovfn = sp.get('cams','fovfn',fallback=None)
#%% manual override flux file
        try:
            self.Jfwdh5 = overrides['Jfwd']
            logging.info('* overriding J0 flux with file: {}'.format(overrides['Jfwd']))
        except (KeyError,TypeError):
            self.Jfwdh5 = None
#%% manual override filter
        try:
            self.opticalfilter = overrides['filter'].lower()
            logging.info('* overriding filter choice with:',overrides['filter'])
        except (KeyError,TypeError):
            self.opticalfilter = sp['transcar']['opticalFilter'].lower()
#%% manual override minimum beam energy
        try:
            self.minbeamev = float(overrides['minev'])
            logging.info('* minimum beam energy set to: {}'.format(overrides['minev']))
        except (KeyError,TypeError): #use spreadsheet
            self.minbeamev = sp.getfloat('transcar','minbeamev',fallback=0)

#%% fit method
        try:
            self.optimfitmeth = overrides['fitm']
            logging.info('* setting fit method to {}'.format(overrides['fitm']))
        except (KeyError,TypeError):
            self.optimfitmeth = sp.get('recon','OptimFluxMethod',fallback='').lower()

        self.optimmaxiter = sp.getint('recon','OptimMaxiter',fallback=None)
#%% force compute ell
        try:
            if overrides and overrides['ell']:
                self.loadfwdL = False
            else:
                self.loadfwdL = True
        except KeyError:
                self.loadfwdL = True
#%% setup plotting
#        self.plots = {}
#
#        if 'optim' in makeplot:
#            self.plots['optim'] = ('bnoise','best','pest','phiest','pest_1d','phiest_1d')
#        if 'fwd' in makeplot:
#            self.plots['fwd'] =   ('bnoise','bfwd','pfwd','phifwd','pfwd_1d','phifwd_1d')
#%%how many synthetic arcs are we using
        self.nArc = len(ap) # number of dict entries
        self.nTimeSlice = ntimeslice
#%% transcar
        self.lambminmax = (splitconf(sp,('sim','lambdamin')),
                           splitconf(sp,('sim','lambdamax'))) #for plotting only

        self.useztranscar = sp.getboolean('transcar','UseTCz',fallback=None)
        self.loadver      = sp.getboolean('transcar','loadVER',fallback=None)
        self.loadverfn    = sp.get('transcar','verfn',fallback=None)
        self.bg3fn =     (self.rootdir/sp['sim']['BG3transFN']).expanduser()
        self.windowfn =  (self.rootdir/sp['sim']['windowFN']).expanduser()
        self.qefn =      (self.rootdir/sp['sim']['emccdQEfn']).expanduser()
        self.transcarev =   (self.rootdir/sp['transcar']['BeamEnergyFN']).expanduser()
        self.transcarutc = parse(sp['transcar']['tReq'])
        self.excratesfn =  sp['transcar']['ExcitationDATfn'] #NO ROOTDIR!
        self.transcarpath = (self.rootdir/sp['sim']['TranscarDataDir']).expanduser()
        self.reactionfn =   (self.rootdir/sp['transcar']['reactionParam']).expanduser()
        self.transcarconfig = sp['transcar']['simconfig']

        self.minflux= sp.getfloat('recon','minflux',fallback=0.)

        self.reacreq = ()
        if sp.getboolean('transcar','metastable',fallback=None): self.reacreq += 'metastable',
        if sp.getboolean('transcar','atomic',fallback=None): self.reacreq += 'atomic',
        if sp.getboolean('transcar','N21NG',fallback=None): self.reacreq += 'n21ng',
        if sp.getboolean('transcar','N2meinel',fallback=None): self.reacreq += 'n2meinel',
        if sp.getboolean('transcar','N22PG',fallback=None): self.reacreq += 'n22pg',
        if sp.getboolean('transcar','N21PG',fallback=None): self.reacreq += 'n21pg',

        self.realdata = sp.getboolean('sim','realdata',fallback=None)
        if self.realdata:
            self.realdatapath = Path(sp['cams']['ActualDataDir']).expanduser()

        self.raymap = sp['cams']['RayAngleMapping'].lower()

        self.downsampleEnergy = sp.getfloat('transcar','downsampleEnergy',fallback=False)

        if self.raymap == 'astrometry':
            logging.info('Using ASTROMETRY-based per-pixel 1D cut mapping to 2D model space')
        elif self.raymap == 'arbitrary':
            logging.info('Using arbitrary linear per-pixel 1D cut mapping to 2D model space')
        else:
            raise ValueError('Unknown Ray Angle Mapping method: {}'.format(self.raymap))

        self.cal1dpath = (self.rootdir/sp['cams']['cal1Ddir']).expanduser()

        try: #override
            if overrides['treq'] is not None: #must be is not None for array case
                if isinstance(overrides['treq'][0],(datetime,float,int)):
                    pass
                elif isinstance(overrides['treq'][0],str):
                    overrides['treq'] = [parse(t).timestamp() for t in overrides['treq']]

                self.treqlist = overrides['treq']

            else:
                self.startutc = parse(sp.get('cams','reqStartUT',fallback=None)) #not fallback=''
                self.stoputc =  parse(sp.get('cams','reqStopUT', fallback=None))
        except (KeyError,TypeError):
            try:
                self.startutc = parse(sp.get('cams','reqStartUT',fallback=None))
                self.stoputc =  parse(sp.get('cams','reqStopUT', fallback=None))
            except (KeyError,ValueError):
                pass


        self.timestepsperexp = sp.getint('sim','timestepsPerExposure')
        #%% recon
        self.artinit =    sp.get('recon','initVector',fallback='').lower()
        self.artmaxiter = sp.getint('recon','maxIter',fallback=0)

        self.artlambda = sp.getfloat('recon','lambda',fallback=None)
        self.artstop =   sp.get('recon','stoprule',fallback=None)
        self.arttau =    sp.getfloat('recon','MDPtauDelta',fallback=None)

    def setupFwdXZ(self,sp):
        Fwd = {}
        self.fwd_xlim = (sp.getfloat('fwd','XminKM'), sp.getfloat('fwd','XmaxKM'))
        self.fwd_zlim = (sp.getfloat('fwd','ZminKM'), sp.getfloat('fwd','ZmaxKM'))
        self.fwd_dxKM = sp.getfloat('fwd','XcellKM')


        if self.useztranscar:
            Fwd['x'] = makexzgrid(self.fwd_xlim, None, self.fwd_dxKM, None)[0]
            zTranscar = getaltgrid(self.rootdir/sp['transcar']['altitudePreload'])
            Fwd['z'] = zTranscar[ (self.fwd_zlim[0] < zTranscar) & (zTranscar < self.fwd_zlim[1]) ]
        else:
            self.fwd_dzKM = sp['fwd']['ZcellKM']
            (Fwd['x'], Fwd['z']) = makexzgrid(self.fwd_xlim, self.fwd_zlim, self.fwd_dxKM, self.fwd_dzKM)
        if Fwd['z'] is None:
            raise ValueError('You must specify zmax zmin zcell when not using transcar altitudes in XLS')

        assert Fwd['x'].ndim == Fwd['z'].ndim ==1
        Fwd['sx'] = Fwd['x'].size #this is a vector
        Fwd['sz'] = Fwd['z'].size #this is a vector

        #maximum number of grid elements in a ray! This number is arrived at 'graphically'
        #by considering the number of cells touched by a ray at a 45 deg. angle at the
        # corner of a square cell block
        Fwd['maxNell'] = ( 2*ceil( hypot( Fwd['sx'], Fwd['sz'] ) ) - 1 ).astype(int)

        if Fwd['sx'] * Fwd['sz'] > 10e6:
            logging.warning('Fwd Grid size seems excessive at more than 10 million cells')

        return Fwd

    def maketind(self,timeInds):
        if self.nTimeSlice == 0:
            raise ValueError('zero frame indices were specified to process')

        if timeInds is None:
            timeInds = arange(self.nTimeSlice) #NOT range!

        timeInds = atleast_1d(timeInds) #necessary for next indexing step

        return timeInds[timeInds<self.nTimeSlice] #(it's <, not <=) slice off commond line requests beyond number of frames

    def getEllHash(self,sp,x,z):
        EllCritParams =  [x, z,
                          fromstring(sp['cam']['nCutPix'],dtype=int,sep=','),
                          fromstring(sp['cam']['FOVmaxLengthKM'],dtype=float,sep=','),
                          sp['cams']['RayAngleMapping'],
                          sp.getfloat('fwd','XcellKM'), sp.getfloat('fwd','XminKM'),
                          sp.getfloat('fwd','XmaxKM'),
                          sp['sim']['EllIs'].lower(),
                          sp.getboolean('transcar','UseTCz')]

        if not self.useztranscar: #FIXME maybe we should always consider these for best safety
            EllCritParams.extend([sp.getfloat('fwd','ZcellKM'),
                                  sp.getfloat('fwd','ZminKM'), sp.getfloat('fwd','ZmaxKM') ])

        if self.raymap == 'arbitrary':
            EllCritParams.extend([fromstring(sp['cam']['boresightElevDeg'],dtype=float,sep=','),
                                  fromstring(sp['cam']['FOVdeg'],dtype=float,sep=',')])

        # get a long string from a mix of numbers and strings
        EllTxt = ''.join([str(n) for n in EllCritParams]) #so fast!
        ellHashed = md5( EllTxt.encode('utf-8') ).hexdigest()

        return ''.join(('Ell_', ellHashed,'.h5'))


def makexzgrid(xLim,zLim,dxKM,dzKM):
    # setup grid
    #it's arange, not range()
    xKM = arange(xLim[0], xLim[1] + dxKM, dxKM, dtype=float) #horizontal sample locations

    if xKM[-1]>xLim[1]:
        xKM = xKM[:-1] #discard last element that was > xLim[1]

    if zLim and isfinite(zLim).all() and dzKM and isfinite(dzKM):
        zKM = arange(zLim[0], zLim[1] + dxKM, dzKM, dtype=float) #vert sample locations
        if zKM[-1]>zLim[1]:
            zKM = zKM[:-1] #discard last element that was > xLim[1]
    else:
        zKM = None
    return xKM,zKM

