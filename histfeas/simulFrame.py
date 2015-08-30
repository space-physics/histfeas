"""
functions for loading HST real camera raw data
 michael hirsch 2014, ported from Matlab code

INPUT FILE FORMAT: intended for use with "DMCdata" raw format, 4-byte
 "footer" containing frame index (must use typecast)
"""
from __future__ import division,absolute_import
from time import time
from numpy import arange, empty, asarray, uint16, rot90, fliplr, flipud
from dateutil import parser
from dateutil.relativedelta import relativedelta
import calendar
from scipy.interpolate import interp1d
# local
import histutils.rawDMCreader as rdr
from .get1Dcut import get1Dcut

def getSimulData(sim,cam,makeplot,progms,verbose=0):
#%% synchronize
    cam,sim = HSTsync(sim,cam,verbose)
#%% load 1-D cut slices into keogram array
    cam,rawdata = HSTframeHandler(sim,cam,makeplot,progms,verbose)
    return cam,rawdata,sim

def HSTsync(sim,cam,verbose):

    try:
        reqStart = parser.parse(sim.startutc)
        reqStop = parser.parse(sim.stoputc)
    except AttributeError: #no specified time
        reqStart = parser.parse("1970-01-01T00:00:00Z") #arbitrary time in the past
        reqStop =  parser.parse("2100-01-01T00:00:00Z")#arbitrary time in the future

#%% get more parameters per used camera
    for c in cam:
        cam[c].ingestcamparam(sim)
#%% determine mutual start/stop frame
# FIXME: assumes that all cameras overlap in time at least a little.
# we will play only over UTC times for which both sites have frames available
    mutualStart = max( [cam[c].startUT for c in cam] ) #who started last
    mutualStop =  min( [cam[c].stopUT  for c in cam] )   # who ended first
#%% make playback time steps
# based on the "simulated" UTC times that do not correspond exactly with either camera, necessarily.
#TODO check for off-by-one
#FIXME use UT1_unix time, will simplify this and following sections
    alltReq = [mutualStart] #makes a list so we can append
    while alltReq[-1] < (mutualStop - relativedelta(seconds=sim.kineticSec)):
        alltReq.append( alltReq[-1] + relativedelta(seconds=sim.kineticSec) )

    nMutRawFrame = len(alltReq)  # NOT alltReq.size--it's a list as it should be!

    print('{} mutual frames available from {} to {}'.format(nMutRawFrame,mutualStart,mutualStop))

#%% adjust start/stop to user request
    alltReqAdj = asarray([t for t in alltReq if t>reqStart and t<reqStop ]) #keep greater than start time
    nMutSim = alltReqAdj.size
    if verbose > 0:
        print('Per user specification, analyzing {} frames from {} to {}'.format(nMutSim,alltReqAdj[0],alltReqAdj[-1]) )
#%% use *nearest neighbor* interpolation to find mutual frames to display.
#   sometimes one camera will have frames repeated, while the other camera
#   might skip some frames altogether
    alltReqUnix = asarray([calendar.timegm(t.utctimetuple()) + t.microsecond / 1e6 for t in alltReqAdj ])

    for c in cam:
        ft = interp1d(cam[c].tCamUnix,
                      arange(cam[c].nFrame,dtype=int),
                      kind='nearest')
        cam[c].pbInd = ft(alltReqUnix).astype(int) #these are the indices for each time (the slower camera will use some frames twice in a row)

    sim.alltReq = alltReqAdj
    sim.nTimeSlice = alltReqAdj.size

    return cam,sim

def HSTframeHandler(sim,cam,makeplot,progms,verbose=0):
#%% load 1D cut coord
    cam = get1Dcut(cam,makeplot,progms,verbose)

#%% use 1D cut coord
    if verbose>0:
        print('frameHandler: Loading and 1-D cutting data...')
    tic = time()
    rawdata = {}
    for c in cam:
        nProcFrame = cam[c].pbInd.size #should be the same for all cameras! FIXME add assert

        keo = empty( ( cam[c].nCutPix, len(cam[c].pbInd) ),dtype=uint16,order='F') #1-D cut data
        tKeo = empty( nProcFrame,dtype=object) #datetime of each frame
        #yes rawdata is order C!
        rawdata[c] = empty( ( nProcFrame, cam[c].SuperX, cam[c].SuperY),
                             dtype=uint16,order='C')
        with open(cam[c].fnStemCam, 'rb') as fid:
            finf = {'bytesperframe':cam[c].BytesPerFrame,
                    'pixelsperimage':cam[c].PixelsPerImage,
                    'nmetadata':cam[c].Nmetadata,
                    'superx':cam[c].SuperX,
                    'supery':cam[c].SuperY
                    }

            for j,iFrm in enumerate(cam[c].pbInd):

                #FIXME compare rawFrameInd with truly requested frame to catch off-by-one errors
                frame,rawFrameInd = rdr.getDMCframe(fid,iFrm,finf,verbose=-1)
                #print(frame.flags)
                #print(iFrm)
                if cam[c].transpose:
                    frame = frame.T
                # rotate -- note if you use origin='lower', rotCCW -> rotCW !
                if cam[c].rotCCW != 0:
                    frame = rot90(frame,k=cam[c].rotCCW)
                # flip
                if cam[c].flipLR:
                    frame = fliplr(frame)
                if cam[c].flipUD:
                    frame = flipud(frame)
                # declare frame UTC time based on raw Index, start time, and kinetic rate
                tKeo[j] = ( cam[c].startUT +
                           relativedelta(seconds= (rawFrameInd - cam[c].firstFrameNum) * cam[c].kineticSec ) +
                           relativedelta(seconds = cam[c].timeShiftSec) )
                #%% do pixel cutting
                keo[:,j] = frame[cam[c].cutrow,cam[c].cutcol]
                #store raw frame for playback synchronized of raw video
                rawdata[c][j,:,:] = frame

        #assign slice & time to class variables
        cam[c].keo = keo
        cam[c].tKeo = tKeo

    if verbose >0: print('DONE  in {:.2f} seconds.'.format(time() - tic))
    return cam,rawdata
