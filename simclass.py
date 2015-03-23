from __future__ import print_function, division
from numpy import asarray,where,arange,isfinite,ceil,hypot
import h5py
from os.path import join
#
import sys
sys.path.append('../transcar-utils')
from readionoinit import getaltgrid

class Sim:

    def __init__(self,sp,cp,ap,overrides,progms,dbglvl):
        self.dbglvl=dbglvl
        #%% how many cameras in use, and which ones?
        usecamreq = asarray(overrides['cam'])
        if usecamreq[0] is not None: #override spreadsheet
            print('* Overriding XLS parameters, using cameras:',usecamreq)
            for ci,civ in enumerate(cp.loc['useCam']): # this might be a silly indexing method, but works
                if (ci==usecamreq).any():
                    cp.ix['useCam',ci] = 1 #do not use boolean, it screws up other rows
                else:
                    cp.ix['useCam',ci] = 0 #do not use boolean, it screws up other rows
        self.useCamBool = cp.loc['useCam'].values.astype(bool)
        self.useCamInd = where(self.useCamBool)[0] # we want first item in tuple

        if usecamreq[0] is not None:
            assert (self.useCamInd == usecamreq).all()

        self.nCam = int(self.useCamBool.sum())

        self.allCamInd = cp.ix['useCam',:].values.astype(bool)
        self.nCutPix = int(cp.ix['nCutPix',0]) #FIXME someday allow different # of pixels..
        self.allCamXkm = cp.ix['Xkm',:].values.astype(float)
        self.allCamZkm = cp.ix['Zkm',:].values.astype(float)

        #%% camera position override
        self.camxreq = overrides['camx']
        if self.camxreq[0] is not None and len(self.camxreq) != self.nCam:
            exit('*** must specify same number of x-loc and used cameras')
        #%% manual override flux file
        if overrides['Jfwd'] is not None:
            print('* overriding J0 flux with file: ' + overrides['Jfwd'])
            self.Jfwdh5 = overrides['Jfwd']
        else:
            self.Jfwdh5 = None
        #%% manual override filter
        if overrides['filter'] is not None:
            print('* overriding filter choice with:',overrides['filter'])
            #sp.loc['OpticalFilter','Transcar'] = overrides['filter']
            self.opticalfilter = overrides['filter'].lower()
        else:
            self.opticalfilter = sp.at['OpticalFilter','Transcar'].lower()
        #%% manual override minimum beam energy
        if overrides['minev'] is not None:
            print('* minimum beam energy set to:',overrides['minev'])
            #sp.loc['minBeameV','Transcar']  = overrides['minev']
            self.minbeamev = overrides['minev']
        else:
            mbe = sp.at['minBeameV','Transcar']
            if isfinite(mbe):
                self.minbeamev = mbe
            else:
                self.minbeamev = 0.
        #%% fit method
        if overrides['fitm'] is not None:
            print('* setting fit method to', overrides['fitm'])
            #sp.loc['OptimFluxMethod','Recon'] = overrides['fitm']
            self.optimfitmeth = overrides['fitm']
        else:
            self.optimfitmeth = sp.at['OptimFluxMethod','Recon']

        self.optimmaxiter = sp.at['OptimMaxiter','Recon']
        #%% force compute ell
        if overrides['ell']:
            #sp.loc['saveEll','Sim'] = 1 #we'll save to out/date directory!
            #sp.loc['loadEll','Sim'] = 0
            self.savefwdL = True
            self.loadfwdL = False
        else:
            self.savefwdL = sp.at['saveEll','Sim']
            self.loadfwdL = sp.at['loadEll','Sim']
        #%%how many synthetic arcs are we using
        if ap is not None:
            self.useArcBool = ap.loc['useArc'].values.astype(bool)
            self.useArcInd = where(self.useArcBool)[0] # we want first item in tuple
            self.nArc = int(self.useArcBool.sum())
        else:
            self.nArc = 0
#%% transcar
        self.useztranscar = sp.at['UseTCz','Transcar'] == 1
        self.loadver = sp.at['loadVER','Transcar'] == 1
        self.loadverfn = sp.at['verfn','Transcar']
        self.bg3fn = sp.at['BG3transFN','Sim']
        self.windowfn = sp.at['windowFN','Sim']
        self.qefn =sp.at['emccdQEfn','Sim']
        self.transcarev = sp.at['BeamEnergyFN','Transcar']
        self.transcarutc = sp.at['tReq','Transcar']
        self.excratesfn = sp.at['ExcitationDATfn','Transcar']
        self.transcarpath = sp.at['TranscarDataDir','Sim']
        self.reactionfn = sp.at['reactionParam','Transcar']
        self.transcarconfig = sp.at['simconfig','Transcar']

        self.reacreq = ()
        if sp.at['metastable','Transcar'] == 1: self.reacreq += 'metastable',
        if sp.at['atomic','Transcar'] == 1: self.reacreq += 'atomic',
        if sp.at['N21NG','Transcar'] == 1: self.reacreq += 'n21ng',
        if sp.at['N2Meinel','Transcar'] == 1: self.reacreq += 'n2meinel',
        if sp.at['N22PG','Transcar'] == 1: self.reacreq += 'n22pg',
        if sp.at['N21PG','Transcar'] == 1: self.reacreq += 'n21pg',

        self.realdata = sp.at['useActualData','Sim'] == 1
        self.realdatapath = sp.at['ActualDataDir','Cams']
        self.raymap = str(sp.at['RayAngleMapping','Obs']).lower()

        if sp.loc['downsampleEnergy','Transcar'] >1:
            self.downsampleEnergy = sp.at['downsampleEnergy','Transcar']
        else:
            self.downsampleEnergy = False

        if overrides['ell'] and progms:
            self.FwdLfn = join(progms,self.getEllHash(sp,cp))
        else:
            self.FwdLfn = join('precompute',self.getEllHash(sp,cp))



        if self.raymap == 'astrometry':
            print('Using ASTROMETRY-based per-pixel 1D cut mapping to 2D model space')
        elif self.raymap == 'arbitrary':
            print('Using arbitrary linear per-pixel 1D cut mapping to 2D model space')
        else:
            exit('*** Unknown Ray Angle Mapping method: ' + str(self.raymap))

        self.cal1dpath = sp.at['cal1Ddir','Cams']

        self.startutc = sp.at['reqStartUT','Obs']
        self.stoputc = sp.at['reqStopUT','Obs']
        # make the simulation time step match that of the fastest camera
        self.kineticSec = 1. / (cp.ix['frameRateHz',self.useCamBool]).max()

        #%% recon
        self.artinit = str(sp.at['initVector','ART']).lower()
        try:
            self.artmaxiter = int(sp.at['maxIter','ART'])
        except ValueError:
            self.artmaxiter = 0
        self.artlambda = sp.at['lambda','ART']
        self.artstop = sp.at['stoprule','ART']
        self.arttau = sp.at['MDPtauDelta','ART']

    def setupFwdXZ(self,sp):
        Fwd = {}
        self.fwd_xlim = (sp.at['XminKM','Fwdf'], sp.at['XmaxKM','Fwdf'])
        self.fwd_zlim = (sp.at['ZminKM','Fwdf'], sp.at['ZmaxKM','Fwdf'])
        self.fwd_dxKM = sp.at['XcellKM','Fwdf']


        if self.useztranscar:
            Fwd['x'] = makexzgrid(self.fwd_xlim, None, self.fwd_dxKM, None)[0]
            zTranscar = getaltgrid(sp.at['altitudePreload','Transcar'])
            Fwd['z'] = zTranscar[(zTranscar <=self.fwd_zlim[1]) ] #TODO clips only top
        else:
            self.fwd_dzKM = sp.at['ZcellKM','Fwdf']
            (Fwd['x'], Fwd['z']) = makexzgrid(self.fwd_xlim, self.fwd_zlim, self.fwd_dxKM, self.fwd_dzKM)
        if Fwd['z'] is None:
            exit('*** sanitycheck: You must specify zmax zmin zcell when not using transcar altitudes in XLS')

        assert Fwd['x'].ndim == Fwd['z'].ndim ==1
        Fwd['sx'] = Fwd['x'].size #this is a vector
        Fwd['sz'] = Fwd['z'].size #this is a vector

        #maximum number of grid elements in a ray! This number is arrived at 'graphically'
        #by considering the number of cells touched by a ray at a 45 deg. angle at the
        # corner of a square cell block
        Fwd['maxNell'] = ( 2*ceil( hypot( Fwd['sx'], Fwd['sz'] ) ) - 1 ).astype(int)

        if Fwd['sx'] * Fwd['sz'] > 10e6:
            print('** sanityCheck: Fwd Grid size seems excessive at more than 10 million cells')
        return Fwd

    def maketind(self,timeInds,nTimeSlice):
        if nTimeSlice == 0:
            exit('*** zero frame indices were specified to process, exiting')

        if timeInds is None:
            timeInds = arange(nTimeSlice) #NOT range!

        return timeInds[timeInds<nTimeSlice] #(it's <, not <=) slice off commond line requests beyond number of frames

    def getEllHash(self,sp,cp):
        from hashlib import md5

        EllCritParams =  [cp.ix['Xkm',:].values, cp.ix['Zkm',:].values,
                          cp.ix['nCutPix',:].values, cp.ix['FOVmaxLengthKM',:].values,
                          sp.at['RayAngleMapping','Obs'].lower(),
                          sp.at['XcellKM','Fwdf'],
                          sp.at['XminKM','Fwdf'], sp.at['XmaxKM','Fwdf'],
                          sp.at['EllIs','Sim'].lower(),
                          sp.at['UseTCz','Transcar'] ]

        if not self.useztranscar:
            EllCritParams.extend([sp.at['ZcellKM','Fwdf'],
                                  sp.at['ZminKM','Fwdf'], sp.at['ZmaxKM','Fwdf'] ])

        if self.raymap == 'arbitrary':
            EllCritParams.extend([cp.ix['boresightElevDeg',:].values,
                                  cp.ix['FOVdeg',:].values])

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

    if zLim is not None and isfinite(zLim).all() and dzKM is not None and isfinite(dzKM):
        zKM = arange(zLim[0], zLim[1] + dxKM, dzKM, dtype=float) #vert sample locations
        if zKM[-1]>zLim[1]:
            zKM = zKM[:-1] #discard last element that was > xLim[1]
    else:
        zKM = None
    return xKM,zKM

