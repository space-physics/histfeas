function [Hra, Targ,HST1, HST2, Pix1,Pix2,HST1RADEC,HST2RADEC] = RunAngles6a(varargin)
clc, close all
%%
% *to get 1D slices:*
% *RunAngles6a(1,0,0,0,1)*
% * to save 1D slices indicies:*
% *RunAngles6a(1,0,1,0,1)*
%
%
% Makes line of sight estimates for auroral feature
% untested, may not work
% Michael Hirsch Apr/May 2013
%
% version 5: load target coordinates from PeakFinder program!
% 
% needs to be updated to use point-n-click of auroral features or automatic detection of features
%
% step 0: input parameters
% step 1: get site coords in ECEF
% step 2: find azimuth angles to auroral feature in xy-plane
% step 3: Estimate line-of-sight intersection via least squares
% step 4: plot results
%
% Tested with Matlab R2013a 64-bit under 64-bit Ubuntu 13.04

%% version checks
% check if we can write KML (requires R2013a or newer)
if verLessThan('matlab','8.0')
    error('This program requires Matlab R2012b or newer to run the geometric transforms')
elseif verLessThan('matlab','8.1')
    warning('This program has limited functionality for Matlab version older than R2013a')
end
%check if we can use Mapping Toolbox 
if license('test','map_toolbox')
    mapTB = true; %#ok<*NASGU> %we have the Mapping Toolbox
else
    mapTB = false;
    warning('Mapping Toolbox required')
end



%% step 0: input parameters

P = length(varargin);
if P>0, Hra.doPlots = varargin{1}; else Hra.doPlots = false; end
if P>1, Hra.doWriteKML = varargin{2}; else Hra.doWriteKML = false; end
if P>2, Hra.writePix = varargin{3}; else Hra.writePix = false; end
if P>3, Hra.doPlane = varargin{4}; else Hra.doPlane = false; end
if P>4, Hra.doMFPonly = varargin{5}; else Hra.doMFPonly = false; end

Hra.comet3Dplot = false;
Hra.plateScale = 'Astrometry';%'MHinterp';%'ASK';

AzEldatapath = ['..',filesep,'AzElMeasured'];

HST1.nxPixel = 512; %number of pixels in x-dimension
HST1.nyPixel = 512; %number of pixels in y-dimension
HST1.nPixel = HST1.nxPixel*HST1.nyPixel;

HST2.nxPixel = 512; %number of pixels in x-dimension
HST2.nyPixel = 512; %number of pixels in y-dimension
HST2.nPixel = HST2.nxPixel*HST2.nyPixel;
%==== target observations ==========
% new for ver 5--load from PeakFinder
if ~Hra.doMFPonly
load('Apr14T854stereoIndicies','stereo1','stereo2','both')

xPix1 = stereo1.col;
xPix2 = stereo2.col;

yPix1 = stereo1.row;
yPix2 = stereo2.row;
else
    eyeball2D = [];
    Targ = [];
end
% using coordinates where upper left is x=1,y=1 
% (remember, ImageJ is 0-indexed from upper left!) 
% (and FITS is 1-indexed from lower left!)
%===================================


d2r = pi/180;
r2d = 180/pi;
%Hra.vers =vers;
Hra.rEarth = referenceEllipsoid('wgs84').MeanRadius; %[m]
Hra.PlotLat = [62,68]; % range of earth latitude to plot [deg]
Hra.PlotLon = [-155, -140]; % range of earth longitude to plot [deg]
Hra.NP = 8; %arbitrary, number of verticies in earth plot
Hra.doPlotEarth = true;


% ============== Lat/Lon from ULN-2550 GPS measurements
% -- Alt from Google Earth (rounded alt to nearest 5m)
HST1.lla = [65.1186367, -147.432975,    500]; %lat[deg],lon[deg],alt[m]
HST2.lla = [65.12657,   -147.496908333, 208]; %lat[deg],lon[deg],alt[m] 
% ====================================================

% ==== B-field parameters, from http://www.ngdc.noaa.gov/geomag-web/?model=igrf#igrfwmm
HST1.Bincl = 77.51; %[deg]
HST1.Bdecl = 19.92; %[deg] east

HST2.Bincl = 77.5; %[deg]
HST2.Bdecl = 19.90; %[deg] east
%=================
HST1.Baz = 180 + HST1.Bdecl; %[deg] % yes, Azimuth to MagZenith is 180 PLUS Bdecl
HST1.Bel = HST1.Bincl; %[deg]

HST2.Baz = 180 + HST2.Bdecl; %[deg] % yes, Azimuth to MagZenith is 180 PLUS Bdecl
HST2.Bel = HST2.Bincl; %[deg]
%% plate scaling
addpath(AzEldatapath); 

switch Hra.plateScale
    case 'Astrometry' %using the Astrometry.net program
  addpath('~/HST/calibration/BIAS_SLOPE')
  addpath('~/HST/matlab/coordXforms')
%--------------------------------------------------------------------------
AstroFN1 = '~/HST/calibration/AVG10_HST1'; 
AstroFN2 = '~/HST/calibration/AVG10_HST2';
CalibrationTimeUTC = datenum([2013 04 14 08 54 00]); %time of Ra/Dec estimation with the image plane calibration data
%--------------------------------------------------------------------------
% get coord in Ra/Dec
[~, HST1RADEC] = readAstrometry(AstroFN1,false,false,true,HST1.nxPixel,HST1.nyPixel);
[~, HST2RADEC] = readAstrometry(AstroFN2,true, false,true,HST2.nxPixel,HST2.nyPixel);

[HST1.azDataDeg, HST1.elDataDeg] = RaDec2AzElKoblick(HST1RADEC.RA,HST1RADEC.Decl,...
                                    HST1.lla(1),HST1.lla(2),CalibrationTimeUTC);

[HST2.azDataDeg, HST2.elDataDeg] = RaDec2AzElKoblick(HST2RADEC.RA,HST2RADEC.Decl,...
                                    HST2.lla(1),HST2.lla(2),CalibrationTimeUTC);

% assuming even number of pixels in both sensor dimensions, get center of
% image plane Az/el
HST1.centAzEl(1) = getCenterValue(HST1.azDataDeg);
HST1.centAzEl(2) = getCenterValue(HST1.elDataDeg);

HST2.centAzEl(1) = getCenterValue(HST2.azDataDeg);
HST2.centAzEl(2) = getCenterValue(HST2.elDataDeg);

% get magnetic zenith RA/Dec
% Note: the azimuth is 180 + Declination, Declination a positive number.
% The magnetic zenith appears slightly west of south to these cameras
% geomagnetic position.
[HST1RADEC.magZen(1), HST1RADEC.magZen(2)] = ...
     AzEl2RaDec( HST1.Baz, HST1.Bel,...
                HST1.lla(1),HST1.lla(2),CalibrationTimeUTC);

[HST2RADEC.magZen(1), HST2RADEC.magZen(2)] = ...
     AzEl2RaDec( HST2.Baz, HST2.Bel,...
                HST2.lla(1),HST2.lla(2),CalibrationTimeUTC);

            
            
    case 'ASK'
        HST1RADEC=[];
% get star pointing estimates from April 14, 2013 data using Hanna
% Dahlgren's ASK starfield pointing program (just load the data)


HST1.azDataDeg = HSTazelRawReader([AzEldatapath,filesep,'HST1az.dat'],HST1.nyPixel,HST1.nxPixel);
HST1.elDataDeg = HSTazelRawReader([AzEldatapath,filesep,'HST1el.dat'],HST1.nyPixel,HST1.nxPixel);

HST2.azDataDeg = HSTazelRawReader([AzEldatapath,filesep,'HST2az.dat'],HST2.nyPixel,HST2.nxPixel);
HST2.elDataDeg = HSTazelRawReader([AzEldatapath,filesep,'HST2el.dat'],HST2.nyPixel,HST2.nxPixel);

% assuming even number of pixels in both sensor dimensions, get center of
% image plane Az/el
HST1.centAzEl(1) = getCenterValue(HST1.azDataDeg);
HST1.centAzEl(2) = getCenterValue(HST1.elDataDeg);

HST2.centAzEl(1) = getCenterValue(HST2.azDataDeg);
HST2.centAzEl(2) = getCenterValue(HST2.elDataDeg);

    case 'MHinterp'
        HST1RADEC=[];
[HST1.Faz,HST1.Fel,HST2.Faz,HST2.Fel,both.starMeas] = PlateScaling();
HST1.centAzEl(1) = HST1.Faz(HST1.nxPixel/2,HST1.nyPixel/2);
HST1.centAzEl(2) = HST1.Fel(HST1.nxPixel/2,HST1.nyPixel/2);

HST2.centAzEl(1) = HST2.Faz(HST2.nxPixel/2,HST2.nyPixel/2);
HST2.centAzEl(2) = HST2.Fel(HST2.nxPixel/2,HST2.nyPixel/2);

% for convenience (at expense of memory) throw up the whole grid of az/el
% to allow reusing older data access methods.
% i.e. HST_.azDataDeg holds the "azimuth" for each camera pixel
% and HST_.elDataDeg holds the "elevation" for each camera pixel
[xg, yg] = meshgrid(1:HST1.nxPixel,1:HST1.nyPixel);
HST1.azDataDeg = HST1.Faz(xg,yg);
HST1.elDataDeg = HST1.Fel(xg,yg);

HST2.azDataDeg = HST2.Faz(xg,yg);
HST2.elDataDeg = HST2.Fel(xg,yg);
end


display(['Center of HST1 (Az,El): (',num2str(HST1.centAzEl(1)),',',...
    num2str(HST1.centAzEl(2)),') [deg]'])
display(['Center of HST2 (Az,El): (',num2str(HST2.centAzEl(1)),',',...
    num2str(HST2.centAzEl(2)),') [deg]'])
%================



%% step 1: get site coordinates in ECEF and Az/El/Range (AER)

%get HST1 ECEF coord.
[HST1.x, HST1.y, HST1.z] = geodetic2ecef(HST1.lla(1)*d2r,...
                                         HST1.lla(2)*d2r,...
                                         HST1.lla(3),...
                                    referenceEllipsoid('wgs84'));

%get HST2 ECEF coord.
[HST2.x, HST2.y, HST2.z] = geodetic2ecef(HST2.lla(1)*d2r,...
                                         HST2.lla(2)*d2r,...
                                         HST2.lla(3),...
                                    referenceEllipsoid('wgs84'));

                                
% get az,el,range vector pointing FROM HST1 TO HST2
[D12.az,D12.el,D12.r] = geodetic2aer(HST2.lla(1),HST2.lla(2),HST2.lla(3),...
                                     HST1.lla(1),HST1.lla(2),HST1.lla(3),...
                                     referenceEllipsoid('wgs84'),'degrees');
                                 
% get east,north,up vector pointing FROM HST1 TO HST2
[D12.e, D12.n, D12.u] = geodetic2enu(HST2.lla(1),HST2.lla(2),HST2.lla(3),...
                                     HST1.lla(1),HST1.lla(2),HST1.lla(3),...
                                     referenceEllipsoid('wgs84'),'degrees');
                                
display(['Distance between HST1 and HST2: ',num2str(D12.r/1e3),' km.'])
display(['Azimuth from HST1 to HST2: ',num2str(D12.az,'%0.2f'),' deg.'])        
display(['Elevation Angle from HST1 to HST2: ',num2str(D12.el,'%0.2f'),' deg.'])
%% step 2


% find 1D projection of magnetic field through image plane
if Hra.doPlane || Hra.doMFPonly
[Pix1,Pix2,HST1toMagZen,HST2toMagZen,HST1toPoints,HST2toPoints] =...
               plotMFP2(HST1,HST2,Hra,HST1RADEC,HST2RADEC); %#ok<ASGLU>
end

if ~Hra.doMFPonly
HST1.MaxAuroralSlantRange = 500e3 / sind(HST1.centAzEl(2));
HST2.MaxAuroralSlantRange= 500e3 / sind(HST2.centAzEl(2)); %[m] for least squares sol'n, maximum slantrange to try and find sol'n
% using ImageJ/Matlab matrix coordinates, where upper left is x=1,y=1
switch Hra.plateScale
    case {'Astrometry','ASK'}
%======= HST1 Observations =======
%azimuth from HST1 to target
Targ.Az1 = HST1.azDataDeg(yPix1,xPix1); %[deg]
%elevation from HST1 to target
Targ.El1 = HST1.elDataDeg(yPix1,xPix1);%[deg]
%======= HST2 Observations =======
%azimuth from HST2 to target
Targ.Az2 = HST2.azDataDeg(yPix2,xPix2); %[deg]
 %elevation from HST2 to target
Targ.El2 = HST2.elDataDeg(yPix2,xPix2); %[deg]
%=================================
    case 'MHinterp'
%======= HST1 Observations =======
%azimuth from HST1 to target
Targ.Az1 = HST1.Faz(xPix1,yPix1); %[deg]
%elevation from HST1 to target
Targ.El1 = HST1.Fel(xPix1,yPix1);%[deg]
%======= HST2 Observations =======
%azimuth from HST2 to target
Targ.Az2 = HST2.Faz(xPix2,yPix2); %[deg]
 %elevation from HST2 to target
Targ.El2 = HST2.Fel(xPix2,yPix2); %[deg]
%=================================    
end

%% step 3: Least Squares estimate of closest point of approach of line of sight to auroral feature

[HST1,HST2,Targ] = LSEtarget(HST1,HST2,Targ);


%% step 4: Plot result, try to make KML if Matlab R2013a
xTickStep = 75;
% ================ plot altitude vs. time
figure
plot(both.tReqUT,Targ.lla(:,3)./1e3)
title(['Target Altitude: ',datestr(both.tReqUT(1),'yyyy-mm-dd')])
ylabel('Altitude [km]')
xlabel('time [UTC]')
set(gca,'xtick',both.tReqUT(1:xTickStep:end),'xgrid','on')
datetick('x','keepticks')
%============== plot 2D position in the mutual frame ======================
% just modify the 3D plot of ECEF points in Targ.xyz?
% first, find the plane with MMSE to points
if Hra.doPlane
[~,~,V] = svd([Targ.x-mean(Targ.x),...
               Targ.y-mean(Targ.y),...
               Targ.z-mean(Targ.z)],0);

%now do plot

% Let's look at the middle of the group of target points, in the direction
% normal to the mutual plane connecting HST1, HST2, and mag. zenith

             

% find mean of target points
Targ.midxyz = [mean(Targ.x),mean(Targ.y),mean(Targ.z)];
[Targ.midlla(1),Targ.midlla(2),Targ.midlla(3)] = ecef2geodetic(...
                    Targ.midxyz(1),Targ.midxyz(2),Targ.midxyz(3),...
                    referenceEllipsoid('wgs84'));
Targ.midlla(1:2) = Targ.midlla(1:2)*r2d;

%draw a line normal to plane through this point (numerical, improvised)
% to get necessary angles, take average of ground level HST1 and HST2
% (recall, by using the mutual plane at all we're assuming the angles are
% constant--lots of error to go around!)
eyeball2D.Az = mean([HST1.decl,HST2.decl]); %NOT plus 180! (we're going north)
eyeball2D.El = 180 -(90 + mean([HST1.incl,HST2.incl])); %going perpendicular from plane
eyeball2D.Range = 10e3; %[m] %how far eyeball is from target centroid (arbitrary)
%we want to be perpendicular to this, a little ways away so that we can see
%the 2D picture
[eyeball2D.xyz(1), eyeball2D.xyz(2), eyeball2D.xyz(3)] = aer2ecef(...
        eyeball2D.Az,eyeball2D.El,eyeball2D.Range,...
        Targ.midlla(1)*d2r,Targ.midlla(2)*d2r,Targ.midlla(3),...
        referenceEllipsoid('wgs84'));
                        


Targ.eyeAxScale = 1/1e3;  %arbitrary
eyeball2D.hcometFig = figure;
eyeball2D.hcometAx = axes('parent',eyeball2D.hcometFig);

%draw mutual plane for reference, using a polygon
ArbitraryIndex = 1300;
MutPlaneVerti.x = [HST1.x,HST2.x,HST1toMagZen.x(ArbitraryIndex)];
MutPlaneVerti.y = [HST1.y,HST2.y,HST1toMagZen.y(ArbitraryIndex)];
MutPlaneVerti.z = [HST1.z,HST2.z,HST1toMagZen.z(ArbitraryIndex)];
fill3(MutPlaneVerti.x*Targ.eyeAxScale,MutPlaneVerti.y*Targ.eyeAxScale,MutPlaneVerti.z*Targ.eyeAxScale,[1 0 0],'facealpha',0.5)

eyeball2D.hcometPks = line( nan,nan,nan,'parent',eyeball2D.hcometAx,'marker','*','markersize',10); 
%make axes be fixed 
set(eyeball2D.hcometAx,'xlim',[-2320000, -2312000]*Targ.eyeAxScale,...
             'ylim',[ -1475000 -1470000]*Targ.eyeAxScale,...
             'zlim',[5830000 5842000]*Targ.eyeAxScale)

%==============
%now place eyeball on line normal to mutual plane going through Targ.midpoint
%(1) place eyeball at distant, perpendicular to mutual plane point
campos(eyeball2D.xyz*Targ.eyeAxScale) % ECEF
%(2) point eyeball at target centroid
camtarget(Targ.midxyz*Targ.eyeAxScale) % ECEF
%==========
if Hra.comet3Dplot
eyeball2D.hcometT = title(eyeball2D.hcometT,'HST i=','fontsize',10);
% passing callback
set(eyeball2D.hcometFig,'keypressFcn',{@updateComet3D,eyeball2D,Targ,both})
setappdata(eyeball2D.hcometFig,'i',1)
else
    set(eyeball2D.hcometPks,'xdata',Targ.x*Targ.eyeAxScale,'ydata',Targ.y*Targ.eyeAxScale,'zdata',Targ.z*Targ.eyeAxScale)
end
zlabel(eyeball2D.hcometAx,'ECEF Z [km]')
ylabel(eyeball2D.hcometAx,'ECEF Y [km]')
xlabel(eyeball2D.hcometAx,'ECEF X [km]')
end
%============== plot lat/long (silly) ===================
%{ 
somewhat useless
figure
plot(both.tReqUT,Targ.lla(:,1))
title(['Target Latitude: ',datestr(both.tReqUT(1),'yyyy-mm-dd')])
ylabel('Latitude [Deg]')
xlabel('time [UTC]')
set(gca,'xtick',both.tReqUT(1:75:end),'xgrid','on')
datetick('x','keepticks')

figure
plot(both.tReqUT,Targ.lla(:,2))
title(['Target Longitude: ',datestr(both.tReqUT(1),'yyyy-mm-dd')])
ylabel('Longitude [Deg]')
xlabel('time [UTC]')
set(gca,'xtick',both.tReqUT(1:xTickStep:end),'xgrid','on')
datetick('x','keepticks')
%}

%{
useless
figure
plot(Targ.lla(:,1),Targ.lla(:,2)), hold on
plot(Targ.lla(1,1),Targ.lla(1,2),'*g')
plot(Targ.lla(end,1),Targ.lla(end,2),'*r')
ylabel('Longitude [Deg]')
xlabel('Latitude [Deg]')
%}
%================================
plotBlob(HST1,HST2,Hra,Targ);
end % if ~Hra.doMFPonly
%% cleanup
if nargout==0, clear, end

end %function

function updateComet3D(hFig,event,h,Targ,both)

i = getappdata(hFig,'i');
switch event.Key
    case {'rightarrow','up'},  i = i+1; 
    case {'leftarrow','down'},   i = i-1; 
    case 'pagedown',    i = i-50;
    case 'pageup',      i=i+50;
    case '', %first run
    otherwise, return
end
%control bounds
i(i<1) = 1;
i(i>both.nMutFrame) = both.nMutFrame;

setappdata(hFig,'i',i)


    set(h.cometPks,'xdata',Targ.x(i)*Targ.eyeAxScale,'ydata',Targ.y(i)*Targ.eyeAxScale,'zdata',Targ.z(i)*Targ.eyeAxScale)

    set(h.cometT,'string',['i=',int2str(i)])


end


function plotBlob(HST1,HST2,Hra,Targ)
%% plots
if Hra.doPlots
[xE, yE, zE] = makeEarth(Hra.PlotLat,Hra.PlotLon,Hra.NP,Hra.rEarth,Hra.doPlotEarth);
hEarth = plotEarth(xE,yE,zE);
hold on

[HST1,HST2] = makeRays(HST1,HST2); 
plotRays(HST1,HST2,hEarth)

plot3(Targ.x,Targ.y,Targ.z,'g*')



%[tmpX,tmpY] = ndgrid
%surf(Hra.HSTplane([1,4,7]),Hra.HSTplane([2,5,8]),Hra.HSTplane([3,6,9]))

end

%% write KML

if Hra.doWriteKML
if verLessThan('matlab','8.1'), warning('No KML without Matlab R2013a or newer')
else %vers >=8.1
    
    Decimation = 1; %<update>
    
kmlwrite('angles.kml',Targ.lla(1:Decimation:end,1),...
                      Targ.lla(1:Decimation:end,2),...
                      Targ.lla(1:Decimation:end,3),...
    'AltitudeMode','RelativeToGround','IconScale',0.5,...
    'Name',Hra.kmlName(1:Decimation:end))
kmlwriteline('ray1.kml',HST1.llaRay(:,1),HST1.llaRay(:,2),HST1.llaRay(:,3),...
     'AltitudeMode','RelativeToGround','Name','HST1')
kmlwriteline('ray2.kml',HST2.llaRay(:,1),HST2.llaRay(:,2),HST2.llaRay(:,3),...
     'AltitudeMode','RelativeToGround','Name','HST2')
end
end

end



