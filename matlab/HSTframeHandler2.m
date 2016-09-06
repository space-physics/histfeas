function [cam1,cam2] = HSTframeHandler2(hFig,event,cam1,cam2,both,writerObj,Pix1,Pix2)
% this version allows control of playback by arrow keys

%% handle index argument passing and key presses
t = getappdata(hFig,'t'); %current time index
if both.manualPeaks
HST1pix = getappdata(hFig,'HST1pix'); %pixels clicked by user for HST1
HST2pix = getappdata(hFig,'HST2pix'); %pixels clicked by user for HST2
end
anyMovingKey = {'rightarrow','uparrow','leftarrow','downarrow','pagedown','pageup'};
%---------
switch event.Key
    case {'rightarrow','uparrow'},  t = t+1; 
    case {'leftarrow','downarrow'},   t = t-1; 
    case 'pagedown',    t = t-1/cam1.kineticSec;
    case 'pageup',      t=t+1/cam1.kineticSec;
    case anyMovingKey, %if exists, update clicked cursor data
        if both.manualPeaks
        set(cam1.hmPk,'xdata',HST1pix(t,1),'ydata',HST1pix(t,2)) %update cursor location
        set(cam2.hmPk,'xdata',HST2pix(t,1),'ydata',HST2pix(t,2)) %update cursor location
        end
    case {'x','q','escape'}, close(hFig), return
    case '',           %t=t+1; %first run, handles t-1 of doFeatureTrack
         if both.fbok, decay440(), end
    case 'space', manualPeaks(cam1,cam2,both,t,hFig)
    otherwise, return
end

%control bounds
t(t<1) = 1;
t(t>both.nMutFrame) = both.nMutFrame;

%------------------
% update Handle variables!
setappdata(hFig,'t',t) %stores the present value of "t"
%=================================

% t is loop iteration index
%% read next HST1/2 frames
[cam1Frame,cam1time,frame1Ind] = getFrameData(cam1,both,t,cam1.rotate,1);
[cam2Frame,cam2time,frame2Ind] = getFrameData(cam2,both,t,cam2.rotate,2);
% ======== put frame-by-frame code after this line
%% oblique xcor code goes here (Hanna)

%% Feature tracking
  featureTrackHST(cam1,cam2,both,t)
%% imgradient code May 2013 (optional)
if both.doFrameGrad
    getFrameGradient(cam1,cam1Frame,cam1time,frame1Ind,1);
    getFrameGradient(cam2,cam2Frame,cam2time,frame2Ind,2);
end
%% do frame subtraction (optional)
 frameSub(cam1Frame,cam2Frame,both)
%% (optional) do 1D data cut through plots
[cam1,cam2]=frameCuts(cam1,cam2,both,t,cam1Frame,cam2Frame,Pix1,Pix2);
%% (optional) do Intensity(elev,time) plots (for jls)
%[cam1,cam2,IntEl1,IntEl2]=intElTime(cam1,cam2,both,cam1Frame,cam2Frame,t,IntEl1,IntEl2);
%% (optional) plot peaks from 1D peak finder program
 plotPeaks(cam1,cam2,both,t)
%% let Matlab plots catch up
   %pause(0.01) %arbitrary, used to let Matlab Java graphics "catch up"
%% video/tiff writing  
   videoHandler(cam1,cam2,both,writerObj,frame1Ind,frame2Ind)  
%% manual writing   
if both.doWriteManualPeaksPNG && strcmp(event.Key,'space')
    if length(cam1.fn)>=16
        stem = cam1.fn(1:16);
    else %shorter
        stem = cam1.fn;
    end
    print(cam1.figH,...
    [both.ManualPeakspngDir,filesep,stem,...
           '_HST1-',int2str(frame1Ind),...
           '_HST2-',int2str(frame2Ind),...
           '.png'],...
           '-dpng','-r0') 
end
end

function manualPeaks(cam1,cam2,both,t,hFig)


        if both.doFrameCuts && both.plotPeaks %use 1D cuts only, with ***PRELOADED*** PEAK DATA 
        xPix1 = cam1.stereo.col(~isnan(cam1.stereo.col(:,t)),t);
        yPix1 = cam1.stereo.row(~isnan(cam1.stereo.row(:,t)),t);
        xPix2 = cam2.stereo.col(~isnan(cam2.stereo.col(:,t)),t);
        yPix2 = cam2.stereo.row(~isnan(cam2.stereo.row(:,t)),t);
        end
        
        if both.manualPeaks % using DataCursors, one per camera
HST1pix = getappdata(hFig,'HST1pix'); %pixels clicked by user for HST1
HST2pix = getappdata(hFig,'HST2pix'); %pixels clicked by user for HST2
Targout = getappdata(hFig,'Targout'); %target data
TargAlt = getappdata(hFig,'TargAlt'); % estimated Target altitude [m]

            axes(cam1.pbAx) %set current axis to HST1
            display('Click in left HST1 pane, or press <Enter> to reuse existing point')
            if both.fbok, singSine(1200), end
            [cam1.newCursX,cam1.newCursY] = ginput(1);
            if ~isempty(cam1.newCursX) && gca == cam1.pbAx %user clicked new point in correct axes
            set(cam1.hmPk,'xdata',cam1.newCursX,'ydata',cam1.newCursY) %update cursor location
            else
                display('Reused existing point:')
            end
            
            xPix1 = fix(get(cam1.hmPk,'xData')); 
            yPix1 = fix(get(cam1.hmPk,'yData'));
            display(['HST1: (xPix1,yPix1) =(',num2str(xPix1),',',num2str(yPix1),')'])
            display('-------------------------------')
            
            display('Click in right HST2 pane, or press <Enter> to reuse existing point')
            if both.fbok, singSine(660), end
            axes(cam2.pbAx) %set current axis to HST2
            [cam2.newCursX,cam2.newCursY] = ginput(1);
            if ~isempty(cam2.newCursX) && gca == cam2.pbAx %user clicked new point in correct axes
            set(cam2.hmPk,'xdata',cam2.newCursX,'ydata',cam2.newCursY) %update cursor location
            else
                display('Reused existing point:')
            end
            
            xPix2 = fix(get(cam2.hmPk,'xData')); 
            yPix2 = fix(get(cam2.hmPk,'yData'));
            display(['HST2: (xPix2,yPix2) =(',num2str(xPix2),',',num2str(yPix2),')'])
            display('-------------------------------')
            
        

% using ImageJ and/or Matlab matrix coordinates, where upper left is x=1,y=1

switch both.plateScale
    case 'ASK'
%======= HST1 Observations =======
%azimuth from HST1 to target
Targ.Az1 = cam1.azDataDeg(yPix1,xPix1); %[deg]
%elevation from HST1 to target
Targ.El1 = cam1.elDataDeg(yPix1,xPix1);%[deg]
%======= HST2 Observations =======
%azimuth from HST2 to target
Targ.Az2 = cam2.azDataDeg(yPix2,xPix2); %[deg]
 %elevation from HST2 to target
Targ.El2 = cam2.elDataDeg(yPix2,xPix2); %[deg]
%=================================
    case 'MHinterp'
%======= HST1 Observations =======
%azimuth from HST1 to target
Targ.Az1 = cam1.Faz(xPix1,yPix1); %[deg]
%elevation from HST1 to target
Targ.El1 = cam1.Fel(xPix1,yPix1);%[deg]
%======= HST2 Observations =======
%azimuth from HST2 to target
Targ.Az2 = cam2.Faz(xPix2,yPix2); %[deg]
 %elevation from HST2 to target
Targ.El2 = cam2.Fel(xPix2,yPix2); %[deg]
%=================================    
end


%estimate target coordinates
[~,~,Targout{t},TargAlt(t)] = LSEtarget2(cam1,cam2,Targ,both,t);



HST1pix(t,:) = [xPix1,yPix1];
HST2pix(t,:) = [xPix2,yPix2];
setappdata(hFig,'HST1pix',HST1pix)
setappdata(hFig,'HST2pix',HST2pix)
setappdata(hFig,'Targout',Targout)
setappdata(hFig,'TargAlt',TargAlt)


 end
        
   


end

function featureTrackHST(cam1,cam2,both,t)
if both.doFeatureTrack
    
   % Step 1: retreive preview frame data too 
   [cam1LastFrame,cam1LastFrame,frame1LastInd] = getFrameData(cam1,both,t-1,-2,1);
   [cam2LastFrame,cam2LastFrame,frame2LastInd] = getFrameData(cam2,both,t-1,-1,2);
   
   % Step 2: detect features of interest in previous frame and current frame
   metricThreshold = 1000;
   blobs1.prev = detectSURFFeatures(cam1LastFrame, 'MetricThreshold', metricThreshold);
   blobs1.curr = detectSURFFeatures(cam1Frame,     'MetricThreshold', metricThreshold);
   
   blobs2.prev = detectSURFFeatures(cam2LastFrame, 'MetricThreshold', metricThreshold);
   blobs2.curr = detectSURFFeatures(cam2Frame,     'MetricThreshold', metricThreshold);
   
   % Step 3: Select Correspondences Between Points Based on SURF Features
  [features1.prev, validBlobs1.prev] = extractFeatures(cam1LastFrame, blobs1.prev);
  [features1.curr, validBlobs1.curr] = extractFeatures(cam1Frame,     blobs1.curr);
  indexPairs1 = matchFeatures(features1.prev, features1.curr, 'Metric', 'SAD', ...
    'MatchThreshold', matchThreshold);
  
  [features2.prev, validBlobs2.prev] = extractFeatures(cam2LastFrame, blobs2.prev);
  [features2.curr, validBlobs2.curr] = extractFeatures(cam2Frame,     blobs2.curr);
  indexPairs2 = matchFeatures(features2.prev, features2.curr, 'Metric', 'SAD', ...
    'MatchThreshold', matchThreshold);  

  % Retrieve locations of matched points for each image
  matchedPoints1.prev = validBlobs1.prev(indexPairs1(:,1),:);
  matchedPoints1.curr = validBlobs1.curr(indexPairs1(:,2),:);
  
  matchedPoints2.prev = validBlobs2.prev(indexPairs2(:,1),:);
  matchedPoints2.curr = validBlobs2.curr(indexPairs2(:,2),:);
  
   % Step 4. Remove Outliers Using Epipolar Constraints
    numTrials = 10000;
    confidence = 100 - 0.01;
  [fMatrix1, epipolarInliers1, status1] = estimateFundamentalMatrix(...
    matchedPoints1.prev, matchedPoints1.curr, 'Method', 'RANSAC', ...
    'NumTrials', numTrials, 'DistanceThreshold', 0.1, ...
    'Confidence', confidence);
  [fMatrix2, epipolarInliers2, status2] = estimateFundamentalMatrix(...
    matchedPoints2.prev, matchedPoints2.curr, 'Method', 'RANSAC', ...
    'NumTrials', numTrials, 'DistanceThreshold', 0.1, ...
    'Confidence', confidence);  

  % If the function fails to find enough inliers or if either epipole is
  % inside the image, the images cannot be rectified. The function will
  % adjust the parameters and start another iteration.
  if status1 ~= 0 || isEpipoleInImage(fMatrix1, size(cam1Frame)) ...
      || isEpipoleInImage(fMatrix1', size(cam1LastFrame))
   % continue;
  end

  if status2 ~= 0 || isEpipoleInImage(fMatrix2, size(cam2Frame)) ...
      || isEpipoleInImage(fMatrix2', size(cam2LastFrame))
   % continue;
  end
  
   
end

end

function featureTrackSURF(currFrame,prevFrame,t)
%% 2) Collect Interest Points from Each Image

MetricThreshold = 1000;

blobs1 = []; blobs2 = [];

%find blob-like features in both images
while length(blobs1)<10 || length(blobs2)<10
    MetricThreshold = MetricThreshold/2;
blobs1 = detectSURFFeatures(prevFrame, 'MetricThreshold', MetricThreshold,...
        'NumOctaves',5);
blobs2 = detectSURFFeatures(currFrame, 'MetricThreshold', MetricThreshold,...
        'NumOctaves',5);
end
end

function plotPeaks(cam1,cam2,both,t)
if both.plotPeaks
Npks1 = length(cam1.stereo.col(:,t));
Npks2 = length(cam2.stereo.col(:,t));
% do HST1
    set(cam1.pPks,'xdata',cam1.stereo.col(:,t),'ydata',cam1.stereo.row(:,t))
    
           set(cam1.pTxt(1),'position',[cam1.stereo.col(1,t),cam1.stereo.row(1,t)])
if Npks1>1, set(cam1.pTxt(2),'position',[cam1.stereo.col(2,t),cam1.stereo.row(2,t)])
 if Npks1>2, set(cam1.pTxt(3),'position',[cam1.stereo.col(3,t),cam1.stereo.row(3,t)])
 else        set(cam1.pTxt(3),'position',[nan,nan])
 end
else       set(cam1.pTxt(2),'position',[nan,nan])
end

% do HST2
    set(cam2.pPks,'xdata',cam2.stereo.col(:,t),'ydata',cam2.stereo.row(:,t))
        
            set(cam2.pTxt(1),'position',[cam2.stereo.col(1,t),cam2.stereo.row(1,t)])
if Npks2>1, set(cam2.pTxt(2),'position',[cam2.stereo.col(2,t),cam2.stereo.row(2,t)])
 if Npks2>2, set(cam2.pTxt(3),'position',[cam2.stereo.col(3,t),cam2.stereo.row(3,t)])
 else        set(cam2.pTxt(3),'position',[nan,nan])
 end
else       set(cam2.pTxt(2),'position',[nan,nan])
end

end %if plotPeaks
end %function


function [cam1,cam2,IntEl1,IntEl2]=intElTime(cam1,cam2,both,cam1Frame,cam2Frame,t,IntEl1,IntEl2)
% for each frame, sort all the data by elevation and plot.
% data will build up column by column, like a keogram.

if both.doIntElTimePlots
    for i = 1:length(cam1.sortedElRow) %<update> to allow different sized camera images
    IntEl1(i,t) = cam1Frame(cam1.sortedElRow(i),cam1.sortedElCol(i));
    IntEl2(i,t) = cam2Frame(cam2.sortedElRow(i),cam2.sortedElCol(i));
    end
    
%update image
set(cam1.hIntElImg,'cdata',IntEl1)
set(cam2.hIntElImg,'cdata',IntEl2)
    
end

end

function [cam1,cam2]=frameCuts(cam1,cam2,both,t,cam1Frame,cam2Frame,Pix1,Pix2)
% cam_cut is a vector of intensity data along the 1D cut
if both.doFrameCuts && ~both.keypressPlayback
    cam1cut1D = zeros(512,1,'uint16');
    cam2cut1D = zeros(512,1,'uint16');
    for i = 1:512 %<update>
cam1cut1D(i) = cam1Frame(Pix1.LLSrow(i),i);
cam2cut1D(i) = cam2Frame(Pix2.LLSrow(i),i);
    end

    
%% make keogram of 1D cuts
% store each time instant in column vector for best efficiency
cam1.keogram(:,t) = cam1cut1D;
cam2.keogram(:,t) = cam2cut1D;
%update keogram display
set(cam1.hKeoImg,'cdata',cam1.keogram')
set(cam2.hKeoImg,'cdata',cam2.keogram')

end

end

function [Gx,Gy,Gmag,Gdir] = getFrameGradient(cam,camFrame,camTime,frameInd,camNum)

%denoise with adaptive Wiener filtering
dnFrame = wiener2(camFrame);

%first get dx,dy
[Gx,Gy] = imgradientxy(dnFrame);
[Gmag,Gdir] = imgradient(Gx,Gy);

    set(cam.gradImgH, 'cdata',Gmag)
    set(cam.gradtH,'string',{['HST',int2str(camNum),' GRADIENT: Time:',...
        datestr(camTime,'yyyy-mm-ddTHH:MM:SS.fff'),'UT',...
        ' frame: ',int2str(frameInd)];...
        })



end %function

function [camFrame,camTime,frameInd] = getFrameData(cam,both,t,rotIn,camNum)
if cam.use
    switch cam.ext
        case '.tif', camFrame = imread(cam.stem,'index',cam.pbInd(t));
        case '.DMCdata'
            %goto the next frame
            fseek(cam.fid,( cam.pbInd(t) -1 )*cam.BytesPerFrame,'bof');
            %ingest the next frame
            camFrame = fread(cam.fid,[cam.SuperY cam.SuperX],'uint16=>uint16',0,'l');  
            %ingest the next frame's metadata
            metadata = fread(cam.fid,cam.nHead16,'uint16=>uint16',0,'l');
            %squish the next frame's metadata to its actual format
            frameInd = typecast([metadata(2) metadata(1)],'uint32');
        case '.mat'
            camFrame = cam.image(:,:,cam.pbInd(t));
            frameInd = cam.pbInd(t);
    end
    
      camFrame = rot90(camFrame,rotIn);    
      if cam.flipud
          camFrame = flipud(camFrame);
      end
      
      if cam.fliplr
          camFrame = fliplr(camFrame);
      end
      
      camTime = cam.startUT +...
          (double(frameInd-cam.firstFrameNum))*cam.kineticSec/both.day2sec +...
          cam.blessTimeShift/both.day2sec;

         if both.playFrames
    %push to display
    set(cam.imgH, 'cdata',camFrame)
            if ~strcmp(cam.ext,'.mat')
            set(cam.titleH,'string',{['HST',int2str(camNum),': Time:',...
                datestr(camTime,'yyyy-mm-ddTHH:MM:SS.fff'),'UT',...
                ' frame: ',int2str(frameInd)];...
                ['file: ',cam.stem]})
            else
                set(cam.titleH,'string',{['HST',int2str(camNum),...
                    'Synthetic Frame: ',int2str(frameInd)];...
                ['file: ',cam.stem]})
            end
         end %if playFrames
end %if

end %function

function frameSub(cam1Frame,cam2Frame,both)
   %normalize data
if both.doFrameSub
   c1fn = double(cam1Frame)./2^16;
   c2fn = double(cam2Frame)./2^14;
   
   subsFrame = c1fn - c2fn;
   %push to display
   set(both.substImgH, 'cdata',subsFrame)
 
    set(both.substtH,'string','( HST1 - HST2 )')
end   

end

function videoHandler(cam1,cam2,both,writerObj,frame1Ind,frame2Ind) %#ok<INUSL>

   if both.doWriteVid
    currVid = getframe(cam1.figH,[0,0,1320,550]);
    writeVideo(writerObj,currVid);
   end
   
   if both.doWriteTiff
    currVid = getframe(cam1.figH,[0,0,1320,550]);
    currVid = rgb2gray(currVid.cdata);
    imwrite(currVid,[both.dataDir,filesep,cam1.fn(1:16),'.tiff'],'tiff',...
        'WriteMode','append','Compression','deflate')
   end
   
   if both.doWritePNG
       
%        imwrite([cam1Frame,zeros(512,1),cam2Frame],[both.dataDir,'/png/',cam1.fn(1:16),...
%            '_HST1-',int2str(cam1.frameInd),...
%            '_HST2-',int2str(cam2.frameInd),...
%            '.png'],'png',...
%            'Comment',['HST1:',int2str(cam1.frameInd),...
%                       ', HST2:',int2str(cam2.frameInd)])

print(cam1.figH,...
    [both.pngDir,filesep,cam1.fn(1:16),...
           '_HST1-',int2str(frame1Ind),...
           '_HST2-',int2str(frame2Ind),...
           '.png'],...
           '-dpng','-r0')

           %'.pgm'],...
           %'-dpgm','-r0')
   end
end

%{
 function I = gray2rgb(I)
I(I>2^16) = 2^16;
I = I./2^16;



I = cat(3,I,I,I); %http://www.mathworks.com/help/images/converting-between-image-types.html

end
%}