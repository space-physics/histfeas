function [cam1,cam2] = HSTcomputeParams(cam1,cam2,both)


if verLessThan('matlab','8.0')
    error('This program requires Matlab R2012b or newer to run the geometric transforms')
elseif verLessThan('matlab','8.1')
    warning('This program has limited functionality for Matlab version older than R2013a')
end

BitsPerByte = 8; % for the x86 CPU era


%% check input
[cam1.dir, cam1.fn, cam1.ext] = fileparts(cam1.stem);
[cam2.dir, cam2.fn, cam2.ext] = fileparts(cam2.stem);

%% get info for sizes
% Ultra (cam1)
switch cam1.ext
    case '.tif' %this case untested for HST
        display('tif untested for HST')
tmp = imfinfo(cam1.stem);
cam1.nFrame = length(tmp);
cam1.tiffinfo = tmp(1);
cam1.nCol = cam1.tiffinfo.Width;
cam1.nRow = cam1.tiffinfo.Height;

    case '.DMCdata'
info = dir(cam1.stem); %get info about the file
if isempty(info), error(['File ',cam1.stem,' was not found.']), end
cam1.fileSizeBytes = info.bytes;

cam1.xPixels=512;
cam1.yPixels=512;
cam1.xBin = 1; %yikes! should read this from XML
cam1.yBin = 1;
cam1.nHead16 = 2; %# of 16-bit header elements 
cam1.nHeadBytes = 2*cam1.nHead16;
cam1.SuperX = cam1.xPixels/cam1.xBin;
cam1.SuperY = cam1.yPixels/cam1.yBin;
cam1.BitsPerPixel = 16; %not what sensor is, but how many bits the camera sends "down the wire"...

cam1.BytesPerFrame = cam1.SuperX * cam1.SuperY * cam1.BitsPerPixel/BitsPerByte + cam1.nHeadBytes;

if cam1.nFrame ~= cam1.fileSizeBytes/cam1.BytesPerFrame
    warning('filename frame numbers not matching contents of file possibly')
end
cam1.nCol = cam1.SuperX;
cam1.nRow = cam1.SuperY;
cam1.fid = nan;

display('Cam1 parameters hard-coded. Should be reading these from XML instead!! ')
%display(cam1)
%-------------------------------------------------------
    case '.mat' %synthetic data file
load(cam1.stem)
cam1.xPixels = ha.nxPixel;
cam1.yPixels = ha.nyPixel;
cam1.xBin = 1;
cam1.yBin = 1;
cam1.SuperX = cam1.xPixels/cam1.xBin;
cam1.SuperY = cam1.yPixels/cam1.yBin;

cam1.nCol = cam1.SuperX;
cam1.nRow = cam1.SuperY;

cam1.image = image1; %not so efficient!

    otherwise, error(['Unrecognized data file extention ',cam1.ext,'  Legal types are: .mat .tif .DMCdata'])
end %switch

% classic (cam2)
switch cam2.ext
    case '.tif'
          display('tif untested for HST')
tmp = imfinfo(cam2.stem{1});
cam2.nFramePerFile = length(tmp);
cam2.nFrame = cam2.nFramePerFile .* length(cam2.stem); %length of total rollover files
cam2.tiffinfo = tmp(1);
cam2.nCol = cam2.tiffinfo.Width;
cam2.nRow = cam2.tiffinfo.Height;

    case '.DMCdata'
info = dir(cam2.stem); %get info about the file
cam2.fileSizeBytes = info.bytes;

cam2.xPixels=512;
cam2.yPixels=512;
cam2.xBin = 1; %yikes! should read this from XML
cam2.yBin = 1;
cam2.nHead16 = 2; %# of 16-bit header elements 
cam2.nHeadBytes = 2*cam2.nHead16;
cam2.SuperX = cam2.xPixels/cam2.xBin;
cam2.SuperY = cam2.yPixels/cam2.yBin;
cam2.BitsPerPixel = 16; %not what sensor is, but how many bits the camera sends "down the wire"...
cam2.BytesPerFrame = cam2.SuperX * cam2.SuperY * cam2.BitsPerPixel/BitsPerByte + cam2.nHeadBytes;

if cam2.nFrame ~= cam2.fileSizeBytes/cam2.BytesPerFrame
    warning('filename frame numbers not matching contents of file possibly')
end
cam2.nCol = cam2.SuperX;
cam2.nRow = cam2.SuperY;
cam2.fid = nan;

display('Cam2 parameters hard-coded. Should be reading these from XML instead!! ')
%display(cam2)

%-------------------------------------------------------
    case '.mat' %synthetic data file
load(cam2.stem)
cam2.xPixels = ha.nxPixel;
cam2.yPixels = ha.nyPixel;
cam2.xBin = 1;
cam2.yBin = 1;
cam2.SuperX = cam2.xPixels/cam2.xBin;
cam2.SuperY = cam2.yPixels/cam2.yBin;

cam2.nCol = cam2.SuperX;
cam2.nRow = cam2.SuperY;

cam2.image = image2; %not so efficient!
    otherwise, error(['Unrecognized data file extention ',cam2.ext,'  Legal types are: .mat .tif .DMCdata'])
end

%% load site coordinates
d2r = pi/180;

% get star pointing estimates 
% from April 14, 2013 data using Hanna
% Dahlgren's ASK starfield pointing program (just load the data)


addpath(both.AzEldatapath);
addpath(both.stereopath);

switch both.plateScale
    case 'ASK' % not used, very inaccurate
        cam1.azDataDeg = HSTazelRawReader([both.AzEldatapath,filesep,'HST1az.dat'],cam1.yPixels,cam1.xPixels);
        cam1.elDataDeg = HSTazelRawReader([both.AzEldatapath,filesep,'HST1el.dat'],cam1.yPixels,cam1.xPixels);

        cam2.azDataDeg = HSTazelRawReader([both.AzEldatapath,filesep,'HST2az.dat'],cam2.yPixels,cam2.xPixels);
        cam2.elDataDeg = HSTazelRawReader([both.AzEldatapath,filesep,'HST2el.dat'],cam2.yPixels,cam2.xPixels);

        % assuming even number of pixels in both sensor dimensions, get center of
        % image plane Az/el
        cam1.centAzEl(1) = getCenterValue(cam1.azDataDeg);
        cam1.centAzEl(2) = getCenterValue(cam1.elDataDeg);

        cam2.centAzEl(1) = getCenterValue(cam2.azDataDeg);
        cam2.centAzEl(2) = getCenterValue(cam2.elDataDeg);

    case 'MHinterp' % reasonably accurate w/o astrometry.net
        [cam1.Faz,cam1.Fel,cam2.Faz,cam2.Fel,both.starMeas] = PlateScaling();
        cam1.centAzEl(1) = cam1.Faz(cam1.yPixels/2,cam1.xPixels/2);
        cam1.centAzEl(2) = cam1.Fel(cam1.yPixels/2,cam1.xPixels/2);

        cam2.centAzEl(1) = cam2.Faz(cam2.yPixels/2,cam2.xPixels/2);
        cam2.centAzEl(2) = cam2.Fel(cam2.yPixels/2,cam2.xPixels/2);
end

both.maxAltEst = 10e3; %[m] maximum altitude to estimate features
%[m] for least squares sol'n, maximum slantrange to find sol'n
cam1.MaxAuroralSlantRange = both.maxAltEst  / sind(cam1.centAzEl(2));
cam2.MaxAuroralSlantRange = both.maxAltEst  / sind(cam2.centAzEl(2)); 

%get HST1 ECEF coord.
[cam1.x, cam1.y, cam1.z] = geodetic2ecef(cam1.lla(1)*d2r,...
                                    cam1.lla(2)*d2r,...
                                    cam1.lla(3),...
                                    referenceEllipsoid('wgs84'));

%get HST2 ECEF coord.
[cam2.x, cam2.y, cam2.z] = geodetic2ecef(cam2.lla(1)*d2r,...
                                    cam2.lla(2)*d2r,...
                                    cam2.lla(3),...
                                    referenceEllipsoid('wgs84'));




end