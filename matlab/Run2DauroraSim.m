% imagerlocations1(arc.HorizOffsetKM,arc.widthKM,writeMovie)
%
%Routine to investigate optimal station distance for different arc widths
%and separation.
%
% initial program created by Hanna Dahlgren, 14 Feb 2012
% functional enhancements by Michael Hirsch, June/July 2013
%
% INPUTS:
% -------
% arc.HorizOffsetKM: B_per-oriented separation between arcs (km)
% arc.widthKM: B_perp-oriented width of auroral arc (default 0.5km)
% arcHorizOffset: B_perp-oriented offset from HST1 site
% 
% Algorithm:
% Step 0) Intake of user parameters
%
% Step 1) sample-by-sample, generate auroral arc (to be replaced by
% chap_gauss_arc.m)???
%   The variable "intpix" is the 1D intensity across the image array that
%   arises from the path integration through the synthetic arcs (arcs
%   modeled with finite Gaussian width and with Chapman function vertically)
%
% Step 2) Plot results
%
function [image1,image2,ha] = Run2DauroraSim(varargin)
clc
%% Step 0) Intake user parameters
P = length(varargin);
if P>0 && ~isempty(varargin{1}), arc.HorizOffsetKM=varargin{1};
else arc.HorizOffsetKM = [0 0.25 .5 1 -2 -4];
end
if P>1 && ~isempty(varargin{2}), arc.widthKM = varargin{2};       else 
    arc.widthKM = 0.1; %[km]
end
if P>2 && ~isempty(varargin{3}), ha.writeMovie = varargin{3};       else ha.writeMovie = false; end

ha.showBeams = false; %show figure of how pixel path integrals traverse auroral region
ha.movieFN = 'beams.avi';
ha.nArc = size(arc.HorizOffsetKM,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Set your parameters here!         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Focal lengths of cameras%in mm
ha.focallengthMM=[50 50];

ha.dist1to2KM=3; %ground distance (km) to second imager

ha.nxPixel = 512; %number of "x-pixels" in camera
ha.nyPixel = 512; %number of "y-pixels" in camera

ha.rowInd = 1; %future 3D primitive approach!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The focal length results in the following FOV
ha.FOVdeg=[atand(8.2/ha.focallengthMM(1)),atand(8.2/ha.focallengthMM(2))];

%degrees per pixel (assuming sqaure sensor)
ha.pixelscale=[ha.FOVdeg(1)/ha.nxPixel, ha.FOVdeg(2)/ha.nxPixel]; 

%HST1 points up-B, and "tilt" is the angle at which HST2 tips toward HST1 
% SUCH THAT the boresight pixels of HST1 and HST2 cross at 110km altitude
ha.boresightCrossingAltKM = 110; %[km]
ha.tilt=atand(ha.dist1to2KM/ha.boresightCrossingAltKM); 

% bin sizes: 
% B_perp: 10m
% B_parallel: 100m 
% 
% grid extents:
% B_perp: -15km to +15km  <--- 30km/10m = 3000 pixels!
% B_parallel:  100..140km      <--- (140-100)km/100m = 400 pixels!
ha.nxGridAurora = 3000;
ha.nyGridAurora = 400;

%% Step 1) Auroral arc generation, point by point iteration
% FIXME this may be replaced with "chap_gauss_arc" perhaps!

[aurora,ha.dmax] = make2Daurora(ha,arc);
%% Step 2) lines of sight (path integrals) through auroral arc
ha = makeFigures(ha,aurora);

[intpix1,intpix2] = computeLOS(ha,aurora);

%% Step 3) plot results
[image1,image2] =finalPlot(ha,intpix1,intpix2,aurora,arc);

if nargout==0, clear, end %eliminate nuisance console printing

end

function [aurora,dmax] = make2Daurora(ha,arc)

aurora=zeros(ha.nxGridAurora,ha.nyGridAurora);
%lim1=1500-arc.widthKM/2*1000;
%lim2=1500+arc.widthKM/2*1000;

%% 1a) Filling the ionosphere with aurora
for z=1:ha.nyGridAurora %for each sample altitude

%Chapman profile
%z=1:400;
H=50; %in 100s meters, so 50*100 m
Z0 = 110; %[km]
zl=(z-Z0)./H;
v=1*exp(0.5*(1-zl-exp(-zl)));

% store the (inter-arc identical) Chapman fcn intensity at each arc x,y bin location
for iArc = 1:ha.nArc
aurora(ha.nxGridAurora/2 ... center of x-samples (where HST1 boresight live)
      + arc.HorizOffsetKM(ha.rowInd,iArc)*100,... recall x-bin-size is 10m
      z...              at this altitude sample bin
      )=v;                %is assigned Chapman fcn brightness "v"
end
%% 1b) Tapering the emissions on the edges of the aurora
dmax=round(arc.widthKM*100/2);
for d1=1:dmax %for each sample horizontal offset bin
    
    value=v-d1^2/dmax^2; % not a Gaussian width taper
    value(value<0) = 0;

    for iArc = 1:ha.nArc
        aurora(ha.nxGridAurora/2 + arc.HorizOffsetKM(ha.rowInd,iArc)*100 + d1, z) = value;
        aurora(ha.nxGridAurora/2 + arc.HorizOffsetKM(ha.rowInd,iArc)*100 - d1, z) = value;
    end
end %for d1
end %for z
%now we have a 2D slice of simulated aurora, let's do something with it
end

function ha = makeFigures(ha,aurora)

if ha.showBeams
ha.testBeamFig = figure;
ha.testBeamAx = axes('parent',ha.testBeamFig);
ha.testBeamImg = imagesc(1:ha.nyGridAurora,1:ha.nxGridAurora,aurora,'parent',ha.testBeamAx); 

hcb = colorbar('peer',ha.testBeamAx);
set(get(hcb,'ylabel'),'string','Arc Intensity [normalized]')

xlabel(ha.testBeamAx,'altitude bins')
ylabel(ha.testBeamAx,'B_\perp bins')
set(ha.testBeamAx,'nextplot','add')

ha.CurrPix1l = text(nan,nan,'1','horizontalalignment','center',...
    'verticalalignment','middle','units','data','color','w','fontsize',12);
ha.CurrPix1r = text(nan,nan,'1','horizontalalignment','center',...
    'verticalalignment','middle','units','data','color','w','fontsize',12);

ha.CurrPix2l = text(nan,nan,'2','horizontalalignment','center',...
    'verticalalignment','middle','units','data','color','k','fontsize',12);
ha.CurrPix2r = text(nan,nan,'2','horizontalalignment','center',...
    'verticalalignment','middle','units','data','color','k','fontsize',12);

ha.CurrLOS1l = plot([nan,nan],[nan,nan],'w--');
ha.CurrLOS1r = plot([nan,nan],[nan,nan],'w--');

ha.CurrLOS2l = plot([nan,nan],[nan,nan],'k--');
ha.CurrLOS2r = plot([nan,nan],[nan,nan],'k--');

ha.tBt = title(ha.testBeamAx,'','interpreter','none');

% focus the view near the boresight intersection
%set(ha.testBeamAx,'add','xlim',[1 ha.nyGridAurora],'ylim',[ha.nxGridAurora/2-200,ha.nxGridAurora/2+200])

end %showBeam
end

function [intpix1,intpix2] = computeLOS(ha,aurora)
%% 2a) For the central (boresight) pixel
%creating a row of pixels

%initialize variables
intpix1=zeros(ha.nyPixel,1);
intpix2=zeros(ha.nyPixel,1);
int1=0;
int2=0;

if ha.showBeams
%startup movie
if ha.writeMovie, 
    videoObj = VideoWriter(ha.movieFN,'Motion JPEG AVI'); 
    videoObj.Quality = 95;
    open(videoObj);
end



    set(ha.CurrLOS1l,'ydata',[ha.nxGridAurora/2 ha.nxGridAurora/2],...
                  'xdata',[1 ha.nyGridAurora])
    set(ha.CurrLOS2l,'ydata',ha.nxGridAurora/2-ha.dist1to2KM*100+ ...
                            [round((1+1000)*10*tand(ha.tilt)),...
                            round((ha.nyGridAurora+1000)*10*tand(ha.tilt))],...
                  'xdata',[1 ha.nyGridAurora])
end

for j=1:ha.nyGridAurora
    %HST1
    xval=0;
    yval=j;
    xIntInd1 = ha.nxGridAurora/2 + xval;

    int1=int1+aurora(xIntInd1,yval);
    %set(ha.CurrPix1l,'position',[yval,xIntInd1])

    
    %HST2 
    %(j+1000)*10 is the number of 10m B_perp bins 
    xval2=round((j+1000)*10*tand(ha.tilt));
    yval2=j;
    xIntInd2 = ha.nxGridAurora/2 ... start in middle of B_perp bins
        - ha.dist1to2KM*100 ... move over by dist1to2KM*100 10m bins
        + xval2; %  come back according to tilt of HST2
    
    int2=int2+aurora(xIntInd2,yval2);
   % set(ha.CurrPix2l,'position',[yval2,xIntInd2])
    
    %pause(0.005) %let java catchup
end
intpix1(ha.nxPixel/2)=int1;
intpix2(ha.nxPixel/2)=int2;

if ha.showBeams
set(ha.tBt,'string',{['white: HST1, pixel ',int2str(ha.nxPixel/2),'=',num2str(int1)];
                     ['black: HST2, pixel ',int2str(ha.nxPixel/2),'=',num2str(int2)]})
if ha.writeMovie
    vFrame = getframe(ha.testBeamFig);
    writeVideo(videoObj,vFrame);    
end
end

%% 2b) compute the left and right sides of boresight
kIndMax = ha.nxPixel/2-1;
for k=1:kIndMax % k starts at the outermost pixel and moves in to the adjacent-to-boresight pixel
    
int1right=0;
int1left=0;
int2right=0;
int2left=0;

if ha.showBeams
LOS1 = [round((1+1000)*10*tand(ha.pixelscale(1)*k)),...
                        round((ha.nyGridAurora+1000)*10*tand(ha.pixelscale(1)*k))];
                    
LOS2l =[round((1+1000)*10*tand(ha.tilt+ha.pixelscale(2)*k)),...
                   round((ha.nyGridAurora+1000)*10*tand(ha.tilt+ha.pixelscale(2)*k))];
LOS2r =[round((1+1000)*10*tand(ha.tilt-ha.pixelscale(2)*k)),...
                   round((ha.nyGridAurora+1000)*10*tand(ha.tilt-ha.pixelscale(2)*k))];

set(ha.CurrLOS1l,'ydata',ha.nxGridAurora/2 + LOS1,...
               'xdata',[1 ha.nyGridAurora])
set(ha.CurrLOS1r,'ydata',ha.nxGridAurora/2 - LOS1,...
               'xdata',[1 ha.nyGridAurora])
           
set(ha.CurrLOS2l,'ydata',ha.nxGridAurora/2 - ha.dist1to2KM*100 + LOS2l,...
               'xdata',[1 ha.nyGridAurora])
set(ha.CurrLOS2r,'ydata',ha.nxGridAurora/2 - ha.dist1to2KM*100 + LOS2r,...
               'xdata',[1 ha.nyGridAurora])     
end
           
for j=1:ha.nyGridAurora
    %HST1
     %(j+1000)*10 is the number of 10m B_perp bins 
    xval1right=round((j+1000)*10*tand(ha.pixelscale(1)*k)); %
    xval1left =round((j+1000)*10*tand(ha.pixelscale(1)*k));
    yval1=j;
    
    a1right=ha.nxGridAurora/2 + xval1right;
    a1right(a1right > ha.nxGridAurora) = ha.nxGridAurora; %clamps maximum to nxGridAurora
        
    a1left=ha.nxGridAurora/2 - xval1left;
    a1left(a1left < 1)=1;
    
    int1right=int1right + aurora(a1right,yval1);
    int1left=int1left+aurora(a1left,yval1);
    %{
    set(ha.CurrPix1r,'position',[yval1,a1right])
    set(ha.CurrPix1l,'position',[yval1,a1left])
    %}
    
    
    %HST2
     %(j+1000)*10 is the number of 10m B_perp bins 
    xval2right=round((j+1000)*10*tand(ha.tilt+ha.pixelscale(2)*k));
    xval2left =round((j+1000)*10*tand(ha.tilt-ha.pixelscale(2)*k));
    yval2=j;
    
    %better to say a2right = ax1; a2left = ax2;
    ax1=ha.nxGridAurora/2 - ha.dist1to2KM*100 + xval2right;
    ax1(ax1 > ha.nxGridAurora) = 3000;
    
    ax2=ha.nxGridAurora/2 - ha.dist1to2KM*100 + xval2left;
    ax2(ax2 < 1)=1;
    
    int2right=int2right+aurora(ax1,yval2);
    int2left=int2left+aurora(ax2,yval2);
    %{
    set(ha.CurrPix2r,'position',[yval2,ax1])
    set(ha.CurrPix2l,'position',[yval2,ax2])
    pause(0.005) %let java catchup
    %}
    
end
intpix1(256-k)=int1left;
intpix1(256+k)=int1right;
intpix2(256-k)=int2left;
intpix2(256+k)=int2right;

if ha.showBeams
 set(ha.tBt,'string',{['white: HST1, pixels ',int2str(ha.nxPixel/2-k),'=',num2str(int1left)...
                       ' & ',int2str(ha.nxPixel/2+k),'=',num2str(int1right)];
                      ['black: HST2, pixels ',int2str(ha.nxPixel/2-k),'=',num2str(int2left)...
                       ' & ',int2str(ha.nxPixel/2+k),'=',num2str(int2right)]})
     
pause(0.005) %let java catchup
if ha.writeMovie
    vFrame = getframe(ha.testBeamFig);
    writeVideo(videoObj,vFrame);    
end
end %showBeams

end %for k

if ha.writeMovie,   close(videoObj);    end

end

function [image1,image2] = finalPlot(ha,intpix1,intpix2,aurora,arc)

%figure(3)
%clf
%plot(intpix1,'+')
%Creating image1
% image1=zeros(512,512);
% for j=1:512
%     image1(j,:)=intpix1;
% end
%Creating image2
% image2=zeros(512,512);
% for j=1:512
%     image2(j,:)=intpix2;
% end

% ahem
image1 = repmat(intpix1(:)',512,1);
image2 = repmat(intpix2(:)',512,1);

figure


%HST1 image plane
subplot(2,2,2)
imagesc(image1)
shading flat
colorbar
title('Image 1, observing in Mag. Zenith')
axis square

%HST2 image plane
subplot(2,2,4)
imagesc(image2)
shading flat
colorbar
title(sprintf('Image 2 (located %s km south of Im1)',num2str(ha.dist1to2KM)))
axis square
subplot(2,2,[1 3])

%HST1
plot([0,0],[0, 200],'k--') %boresight line
hold on
plot([0 ha.nyGridAurora/2*tand(ha.nxPixel/2*ha.pixelscale(1))],...
     [0 ha.nyGridAurora/2*tand(ha.nxPixel/2*ha.pixelscale(1))/tand(ha.nxPixel/2*ha.pixelscale(1))],'k');

plot([0 ha.nyGridAurora/2*tand(-ha.nxPixel/2*ha.pixelscale(1))],...
     [0 ha.nyGridAurora/2*tand(-ha.nxPixel/2*ha.pixelscale(1))/tand(-ha.nxPixel/2*ha.pixelscale(1))],'k');
%-------------------------
%HST2
%boresight line
plot([-ha.dist1to2KM,...
      ha.nyGridAurora/2*tand(ha.tilt)-ha.dist1to2KM],...
      [0,...
      ha.nyGridAurora/2*tand(ha.tilt)/tand(ha.tilt)],'b--') 

plot([-ha.dist1to2KM,...
      ha.nyGridAurora/2*tand(ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)-ha.dist1to2KM],...
     [0,...
      ha.nyGridAurora/2*tand(ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)/tand(ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)],'b');

plot([-ha.dist1to2KM ha.nyGridAurora/2*tand(-ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)-ha.dist1to2KM],...
     [0 ha.nyGridAurora/2*tand(-ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)/tand(-ha.nxPixel/2*ha.pixelscale(2)+ha.tilt)],'b');

%Plotting the arcs
for iArc = 1:ha.nArc
    x1=ha.nxGridAurora/2 + arc.HorizOffsetKM(ha.rowInd,iArc)*100 - ha.dmax;
    x2=ha.nxGridAurora/2 + arc.HorizOffsetKM(ha.rowInd,iArc)*100 + ha.dmax;
    
    imagesc([arc.HorizOffsetKM(ha.rowInd,iArc) - arc.widthKM/2,...
             arc.HorizOffsetKM(ha.rowInd,iArc) + arc.widthKM/2],...
             [100 140],...
             rot90(aurora(x1:x2,1:ha.nyGridAurora),3)) 
end


axis([-15 15 0 140])
xlabel('Ground distance to north (km)')
ylabel('Altitude (km)')
end
