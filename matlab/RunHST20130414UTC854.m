% a nice auroral event
%
% It's best to copy the raw data files to your PC--the network drive is 
% slow. The program expects the data files to live under:
% WINDOWS:   c:\HSTdata\DataField\2013-04-14\
% Mac/Linux: ~/HSTdata/DataField/2013-04-14/
%
% Press F5 to run this program -- then look at the Figure and press Space
% Bar to allow clicking a point in the left and right panes for
% triangulation!
%-----------------------------
% clicking output
% after you've clicked several point pairs, you can plot altitude vs time
% by pasting into command window:
% TargAlt = getappdata(gcf,'TargAlt'); %target data
% figure,plot(TargAlt/1e3,'r.','markersize',12),ylabel('Altitude [km]'),xlabel('frame #')
%
%
%-------------------------
%
%RunSimulPlayHST -> simulFrameHST -> HSTframeHandler
%
%Cam1 is HST1 (ixon Ultra)
%Cam2 is HST2 (ixon classic)
%
% Michael Hirsch Mar-May 2013
% tested with Matlab R2013a 64-bit Linux
%
% Recommended usage (so that you can manipulate already-loaded data)
%  [cam1,cam2,both,Pix1,Pix2]=RunHST20130414UTC855() 
% 
function [cam1,cam2,both,Pix1,Pix2]=RunHST20130414UTC854(varargin)

P = length(varargin);
if P>0, dataRoot = varargin{1}; 
elseif ispc
    dataRoot = 'C:';  
else %isunix | ismac
    dataRoot = '~';
end

dataRoot = [dataRoot,'/data/2013-04-14/hst/'];

camfn = {'2013-04-14T07-00-CamSer7196_frames_363000-1-369200.DMCdata',...
         '2013-04-14T07-00-CamSer1387_frames_205111-1-208621.DMCdata'};

[cam1, cam2,both] = HSTparams(dataRoot,camfn);

%override
both.reqStartUT = datenum([2013 4 14 8 54 21]); 
%both.reqStopUT = datenum([2013 4 14 8 54 40]); 
% both.reqStartUT = datenum([2013 4 14 8 54 22]); 
 both.reqStopUT = datenum([2013 4 14 8 54 30]); 


[cam1,cam2,both,Pix1,Pix2]=simulFrameHST3(cam1,cam2,both);

%suppress nuisance text if no output arguments

if ~nargout
   clear 
end

end

function [cam1, cam2,both] = HSTparams(dataRoot,camfn)
%% parameters that must (for now) be manually set by the user

both.dataDir = dataRoot;
%-----------------------------------------------------------
%% HST1 user-set parameters
cam1.use = true;
cam1.fn = camfn{1};
cam1.kineticSec = 1/53; % get this from XML file
cam1.fullFileStart = datenum([2013 04 14 06 59 55]); %last line of NMEA file
cam1.firstFrameNum = str2double(regexp(cam1.fn,'(?<=frames_)\d*(?=-1-\d*.DMCdata)','match')); %see filename
cam1.lastFrameNum = str2double(regexp(cam1.fn,'(?<=frames_\d*-1-)\d*(?=.DMCdata)','match')); %see filename
cam1.blessTimeShift = -0.188679245283019; %[seconds]  %compensates for system timing errors (ouch!)
cam1.rotate = 0;
cam1.fliplr = true;
cam1.flipud = false;

%% HST2 user-set parameters
cam2.use = true;
cam2.fn = camfn{2};
cam2.kineticSec = 1/30; %get this from XML file
cam2.fullFileStart = datenum([2013 04 14 07 00 07]); %last line of NMEA file
cam2.firstFrameNum = str2double(regexp(cam2.fn,'(?<=frames_)\d*(?=-1-\d*.DMCdata)','match')); %see filename
cam2.lastFrameNum = str2double(regexp(cam2.fn,'(?<=frames_\d*-1-)\d*(?=.DMCdata)','match')); %see filename
cam2.blessTimeShift = 0; %[seconds] %compensates for system timing errors (ouch!)
cam2.rotate = -1;
cam2.fliplr = false;
cam2.flipud = true;
%-----------------------------------------------------
%% now the rest of the parameters
both.keypressPlayback = true; % <beta> lets user control frame advance
both.playFrames = true; %playback data like movie
%-----------------------
both.doFeatureTrack = false; %not yet working
%---------------------
both.manualPeaks = true; % you manually click the pairs of points--cool!
both.plotClickLive = false; % make windows showing LLA results of target as you click
both.doWriteManualPeaksPNG = false; %make a PNG showing where you clicked superimposed on the image
both.ManualPeakspngDir = 'out';
%---------------------------
% *** keogram, also set "both.keypressPlayback = false, %both.playFrame=true"
both.doFrameCuts = true; %take 1D cuts and plot separately (for stack plot and keogram)
both.plotPeaks = false; %uses 'stereoIndicies.mat' output of HST1DpeakFinder.m
%-----------------------
both.doFrameSub =false; %show cam1-cam2 frames (not really useful)
%-------------------------
both.doWriteVid = false;
both.doWriteTiff = false;
both.doWritePNG = false;
%---------------------------
both.doFrameGrad = false; % take gradient of images
%------------------------------------------------------
both.doIntElTimePlots = false; %plot Intensity(elev,time) for JLS <--takes a long time, may not be useful
%-------------------------------------------
both.diagLSE = false; %extra computations for debugging LSE intersection estimate, KML 
%---------------------------------------------------------
both.plateScale = 'MHinterp';%'ASK'


both.day2sec = 86400;
both.simKineticSec = cam1.kineticSec; %sets the time step of requested playback--suggest setting this equal to the fastest frame rate in use
both.AzEldatapath = ['..',filesep,'AzElMeasured'];
both.stereopath = ['precompute'];
both.plotPeaksFN = 'Apr14T854_stereoIndicies.mat';


% ============== Lat/Lon from ULN-2550 GPS measurements
% -- Alt from Google Earth and SRTM data
cam1.lla = [65.1186367,-147.432975,500]; %lat[deg],lon[deg],alt[m]
cam2.lla = [65.12657,-147.496908333,208]; %lat[deg],lon[deg],alt[m] 
% ====================================================

%-----------------
%these will be overridden by "Run" master program, they're only here in case not set there.
both.reqStartUT = nan; 
both.reqStopUT = nan;
%---------------------
%rotate HST1 data CCW
cam1.rot90ccw = true;

%rotate HST2 data  CCW
cam2.rot90ccw = true;

% data number scaling
if cam1.use
    % data number scaling (arbitrary, to get the best look on your 8-bit PC
    % display)
cam1.minVal = 0; 
cam1.maxVal = 42000;
cam1.minGrad = 0;
cam1.maxGrad = 20e3;

cam1.stem = [both.dataDir,cam1.fn];
cam1.nFrame = cam1.lastFrameNum - cam1.firstFrameNum + 1; 
cam1.startUT = cam1.fullFileStart + ...
               (cam1.firstFrameNum-1)*cam1.kineticSec/both.day2sec;
cam1.stopUT = cam1.fullFileStart +...
                (cam1.lastFrameNum-1)*cam1.kineticSec/both.day2sec;
else
    cam1.startUT = nan; cam1.stopUT = nan; cam1.nFrame = nan;  cam1.stem = '';
end

% HST2 computations
if cam2.use           
cam2.minVal = 950; 
cam2.maxVal = 4000;
cam2.minGrad = 0;
cam2.maxGrad = 1e3;

cam2.stem = [both.dataDir,cam2.fn];
cam2.startUT = cam2.fullFileStart + ...
               (cam2.firstFrameNum-1)*cam2.kineticSec/both.day2sec;
cam2.stopUT = cam2.fullFileStart +...
                (cam2.lastFrameNum-1)*cam2.kineticSec/both.day2sec;
cam2.nFrame = cam2.lastFrameNum - cam2.firstFrameNum + 1;
else
    cam2.startUT = nan; cam2.stopUT = nan; cam2.nFrame = nan; cam2.stem = '';
end

end