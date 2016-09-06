function [cam1, cam2,both,Pix1,Pix2] = simulFrameHST3(cam1,cam2,both)
%
% version 3: allows for loading of synthetic "mat" files from
% HST/Simulation/Run2Daurorasim.m
% 
% Note: simulFrameHST0 -- cam1: HST1 
% cam2: HST2
%
% INPUTS:
% -------
% cam1.stem = filename of camera file(s)
% cam1.kineticSec = original camera kinetic time [seconds]
% cam1.startUT = DATENUM format of camera start time in UT
%
% cam1.stem = filename of camera file(s)
% cam1.kineticSec = original camera kinetic time [seconds]
% cam1.startUT = DATENUM format of camera start time in UT
%
% both.simKineticSec  = kinetic rate of 'simulation time' (secoonds) != 1/both.pbFPS
%
% this function plays back CCD and sCMOS data at the same time, using the
% "nearest neighbor" frame for each--frames can be repeated during playback, 
% especially if camera frame rates are different from each other.
% e.g. if one camera fps is 1/2 the other camera fps, the first camera
% playback will show the same frame twice for every one of the second
% camera
%
% tested with Matlab R2013a 64-bit on Ubuntu 13.04 64-bit 
%
% INPUT FILE FORMAT: intended for use with 16-BIT multipage TIFF grayscale
% OR "DMCdata" raw format
%
% We assume the PC doesn't have enough RAM to load all the files at once,
% so we load and play frame-by-frame
%
% Michael Hirsch Dec 2012 (ver1) / July 2013 (ver3)

% example:
% ccd.stem = '/media/small/BigData/Imaging/2012-12-20/CCD/2012-12-20T22-48-34_frames_317600-1-339300.tif';
% ccd.kineticSec = 0.02711; ccd.startUT = datenum([2012 12 21 1 12 04]);
% cmos.stem = '/media/small/BigData/Imaging/2012-12-20/Neo/X109-8bitogv.tif';
% cmos.kineticSec = 0.03030303; cmos.startUT =datenum([2012 12 21 1 17 43]);
%
% both.simKineticSec = 0.1; %skips frames to speedup playback
% both.reqStartUT = nan; 
% both.reqStopUT = nan;
% [ccd2,cmos2,both2]=simulPlay(ccd,cmos,both)


fclose('all');

%% synchronize
[cam1,cam2,both] = HSTsync3(cam1,cam2,both);
%% setup figures
% we do this after HSTsync fcn since we use frame numbers on plot axes
[cam1,cam2,both,Pix1,Pix2]=makeFigures(cam1,cam2,both);

%% (optional) load az/el tables
[cam1,cam2,IntEl1,IntEl2]=loadAzElPlot(cam1,cam2,both);

%% load 1D peaks (optional)
if both.plotPeaks
   load(both.plotPeaksFN,'stereo1','stereo2')
   cam1.stereo = stereo1;
   cam2.stereo = stereo2;
end
%% playback

%open original raw data files
if strcmp(cam1.ext,'.DMCdata'), cam1.fid = fopen(cam1.stem); end
if strcmp(cam2.ext,'.DMCdata'), cam2.fid = fopen(cam2.stem); end

 
%setup video output
if both.doWriteVid
    writerObj = VideoWriter([both.dataDir,filesep,'myVid.avi'],'Uncompressed AVI');
    open(writerObj);
else
    writerObj = [];
end

%% enable keypress responsiveness
set(cam1.figH,'keypressFcn',{@HSTframeHandler2,cam1,cam2,both,writerObj,Pix1,Pix2})
setappdata(cam1.figH,'t',1) %initializes passed variable 't' to 1
setappdata(cam1.figH,'Targout',{[]}) %initializes passed variable 't' to 1
setappdata(cam1.figH,'HST1pix',nan(both.nMutFrame,2)) %initializes passed variable 'HST1pix'
setappdata(cam1.figH,'HST2pix',nan(both.nMutFrame,2)) %initializes passed variable 'HST2pix'
setappdata(cam1.figH,'TargAlt',nan(both.nMutFrame,1)) %init target altitude [m]
%% decide if going to play all at once, or let user plot frame-by-frame
if ~both.keypressPlayback
    event.Key = 'rightarrow';


    %% setup PNG directory
    if both.doWritePNG
        both.pngDir = [both.dataDir,filesep,'png180ms'];
        mkdir(both.pngDir)
    end

    %{
    %===============
    %temp variables for parfor
    cam1imgH = cam1.imgH;
    cam2imgH = cam2.imgH;
    cam1titleH = cam1.titleH;
    cam2titleH = cam2.titleH;
    cam1tUT = cam1.tUT;
    cam2tUT = cam2.tUT;
    cam1stem = cam1.stem;
    cam2stem = cam2.stem;
    cam1pbInd = cam1.pbInd;
    cam2pbInd = cam2.pbInd;
    cam1BytesPerFrame = cam1.BytesPerFrame;
    cam2BytesPerFrame = cam2.BytesPerFrame;
    cam1SuperY = cam1.SuperY;
    cam2SuperY = cam2.SuperY;
    cam1SuperX = cam1.SuperX;
    cam2SuperX = cam2.SuperX;
    cam1nHead16 = cam1.nHead16;
    cam2nHead16 = cam2.nHead16;
    cam1figH = cam1.figH;
    bothdataDir = both.dataDir;
    cam1fn = cam1.fn;
    %==========
    %}
    hW = waitbar(0);

    %% frame-by-frame processing
    for t = 1:both.nMutFrame

     [cam1, cam2,IntEl1,IntEl2] = HSTframeHandler(cam1,cam2,both,writerObj,t,Pix1,Pix2,IntEl1,IntEl2);
    %{

       %update HST1 frame
      cam1fid = fopen(cam1stem); 
                %goto the next frame
                fseek(cam1fid,( cam1pbInd(t) -1 )*cam1BytesPerFrame,'bof');

                %ingest the next frame
                cam1Frame = fread(cam1fid,[cam1SuperY cam1SuperX],'uint16=>uint16',0,'l');  
                %cam1Frame = rot90(fread(fid1,[cam1.SuperY cam1.SuperX],'uint16=>uint16',0,'l').',1);  
                %ingest the next frame's metadata
                metadata = fread(cam1fid,cam1nHead16,'uint16=>uint16',0,'l');
                %squish the next frame's metadata to its actual format
                cam1frameInd = typecast([metadata(2) metadata(1)],'uint32');
     fclose(cam1fid);

    %      if cam1.rot90ccw
             cam1Frame = rot90(cam1Frame,-2); %rotate 180 degrees
    %      else %unrotated
    %      end

        %push to display
        set(cam1imgH, 'cdata',cam1Frame)
        set(cam1titleH,'string',{['HST1: Time:',...
            datestr(cam1tUT(cam1pbInd(t)),'yyyy-mm-ddTHH:MM:SS.fff'),'UT',...
            ' frame: ',int2str(cam1frameInd)];...
            ['file: ',cam1stem]})

     %======================   
        %update HST2 frams

    cam2fid = fopen(cam2stem); 
                %goto the next frame
                fseek(cam2fid,( cam2pbInd(t) -1 )*cam2BytesPerFrame,'bof');
                 %ingest the next frame
                cam2Frame = fread(cam2fid,[cam2SuperY cam2SuperX],'uint16=>uint16',0,'l');
                %cam2Frame = rot90(fread(fid2,[cam2.SuperY cam2.SuperX],'uint16=>uint16',0,'l').',1);  
                metadata = fread(cam2fid,cam2nHead16,'uint16=>uint16',0,'l');
                cam2frameInd = typecast([metadata(2) metadata(1)],'uint32');
    fclose(cam2fid);

    %     if cam2.rot90ccw

             cam2Frame = rot90(cam2Frame,-1); %rotate 90 deg. CW
    %     else %unrotated        
    %     end

            %push cam2 to display
        set(cam2imgH, 'cdata',cam2Frame)
        %set(cam2.imgH, 'cdata',cam2tframe)
        set(cam2titleH,'string',{['HST2: Time:',...
            datestr(cam2tUT(cam2pbInd(t)),'yyyy-mm-ddTHH:MM:SS.fff'),'UT',...
            ' frame: ',int2str(cam2frameInd)];...
            ['file: ',cam2stem]})
    %% video/tiff writing    



    %        imwrite([cam1Frame,zeros(512,1),cam2Frame],[both.dataDir,'/png/',cam1.fn(1:16),...
    %            '_HST1-',int2str(cam1.frameInd),...
    %            '_HST2-',int2str(cam2.frameInd),...
    %            '.png'],'png',...
    %            'Comment',['HST1:',int2str(cam1.frameInd),...
    %                       ', HST2:',int2str(cam2.frameInd)])
    print(cam1figH,...
        [bothdataDir,'/png/',cam1fn(1:16),...
               '_HST1-',int2str(cam1frameInd),...
               '_HST2-',int2str(cam2frameInd),...
               '.png'],...
               '-dpng','-r0')

    %}

     if ~mod(t,10)
      waitbar(t/both.nMutFrame,hW,[num2str(t/both.nMutFrame*100,'%0.1f'),'%'])
     end
    end %for t

    %% make stack plots

    if both.doFrameCuts
    %x1labels = repmat(cam1.firstFrameNum+cam1.pbInd-1,[512,1])';
    cam1.hfStack = figure('Name','HST1 Stack Plot');
    stackedplot(cam1.keogram,4,20,'x');
    set(cam1.hfStack,'pos',[100,100,560,1050])
    title(['HST1 Stack Plot', datestr(both.tReqUT(1)),'UT to ',...
        datestr(both.tReqUT(end)),'UT'])
    view(-90,75)
    ylabel('1D cut pixel #')
    xlabel('Time')
    %zlabel('Data Numbers')


    cam2.hfStack = figure('Name','HST2 Stack Plot');
    stackedplot(cam2.keogram,4,20,'x');
    set(cam2.hfStack,'pos',[700,100,560,1050])
    view(-90,75)
    title(['HST2 Stack Plot', datestr(both.tReqUT(1)),'UT to ',...
        datestr(both.tReqUT(end)),'UT'])
    ylabel('1D cut pixel #')
    xlabel('Time')
    %zlabel('Data Numbers')

    end %if doCuts
    %% cleanup/shutdown
    try close(hW), end

    if both.doWriteVid
            close(writerObj); 
    end
    %% close files
     if strcmp(cam1.ext,'.DMCdata'), fclose(cam1.fid); end
     if strcmp(cam2.ext,'.DMCdata'), fclose(cam2.fid); end

else %user presses key to advance
    %do first plot
    event.Key = '';
    disp('Press <left arrow> to go back in time, press <right arrow> to go forward in time')
    disp('Use <page down> and <page up> to take correspondingly larger jumps backward/forward')
    disp('----------------')
    disp('<space bar> outputs current triangulation result in Command Window')

    if both.manualPeaks
       disp('****************')
       disp(' You must put a data cursor in each frame for triangulation. ')
       disp('****************')
    end

     HSTframeHandler2(cam1.figH,event,cam1,cam2,both,writerObj,Pix1,Pix2);
end %if
end %function



function [cam1,cam2,IntEl1,IntEl2]=loadAzElPlot(cam1,cam2,both)

if both.doIntElTimePlots
%we'll use these for plot axes, and for retreiving intensity data in
%frameHandler
[cam1.sortedEl,cam1.sortedElRow,cam1.sortedElCol] = sortall(cam1.elDataDeg);
[cam2.sortedEl,cam2.sortedElRow,cam2.sortedElCol] = sortall(cam2.elDataDeg);

%preallocate very-large matricies--this code needs to be improved <update>
IntEl1 = zeros(cam1.xPixels*cam1.yPixels,length(cam1.pbInd),'uint16');
IntEl2 = zeros(cam2.xPixels*cam2.yPixels,length(cam2.pbInd),'uint16');

%update plot axes
% set(cam1.hIntElAx,'ylim',[cam1.sortedEl(1),cam1.sortedEl(end)])
% set(cam2.hIntElAx,'ylim',[cam2.sortedEl(1),cam2.sortedEl(end)])

 %% setup separate keogram figure
if both.doIntElTimePlots
    %overall keogram figure
   both.hfIntEl = figure('Name','Intensity Keogram','pos',[20 20 1200 600]);
   
   %HST1 panel
   cam1.hpIntEl = uipanel('parent',both.hfIntEl,'units','pixels',...
       'pos',[1 1 600 500]);
   cam1.hIntElAx = axes('parent',cam1.hpIntEl);
   cam1.hIntElImg = imagesc(1:both.nMutFrame,cam1.sortedEl,...
                          zeros(cam1.xPixels*cam1.yPixels,both.nMutFrame,'uint16'));
   set(cam1.hIntElAx,'ydir','normal')
   xlabel('t=')
   ylabel('Elevation Angle [deg]')
   title(['HST1 Keogram of elevation ', datestr(both.tReqUT(1)),'UT to ',...
    datestr(both.tReqUT(end)),'UT'])

   cam1.hIntElCB = colorbar('peer',cam1.hIntElAx);
   set(get(cam1.hIntElCB,'ylabel'),'string','Data Numbers')
   
   %fix annoying exponential axis labels
   set(cam1.hIntElAx,'YTickLabel',sprintf('%7d|',get(cam1.hIntElAx,'YTick')))
   
   %HST2 panel
   cam2.hpIntEl = uipanel('parent',both.hfIntEl,'units','pixels',...
       'pos',[610 1 600 500]);
   cam2.hIntElAx = axes('parent',cam2.hpIntEl);
   cam2.hIntElImg = imagesc(1:both.nMutFrame,cam1.sortedEl,...
                          zeros(cam2.xPixels*cam2.yPixels,both.nMutFrame,'uint16'));
   set(cam2.hIntElAx,'ydir','normal')
   title(['HST2 Keogram of elevation ', datestr(both.tReqUT(1)),'UT to ',...
    datestr(both.tReqUT(end)),'UT'])
   xlabel('t=')
   ylabel('Elevation Angle [deg]')
   cam2.hIntElCB = colorbar('peer',cam2.hIntElAx);
   set(get(cam2.hIntElCB,'ylabel'),'string','Data Numbers')
   
   %fix annoying exponential axis labels
   set(cam2.hIntElAx,'YTickLabel',sprintf('%7d|',get(cam2.hIntElAx,'YTick')))
   
end

else
    IntEl1 = [];
    IntEl2 = [];

end

end

function resize_fcn(currentPanelHandle)
%resize axes to fill panel, leaving space for colorbar
set(currentPanelHandle,'units','pixels');
set(gca,'units','pixels');
w_pos = get(currentPanelHandle, 'position');
set(gca, 'position', [2 5 w_pos(3)-100 w_pos(4)-50]);
end

function [cam1,cam2,both,Pix1,Pix2]= makeFigures(cam1,cam2,both)
both.feedbPath = ['..',filesep,'feedback'];
testFB = dir(both.feedbPath);
if ~isempty(testFB)
addpath(both.feedbPath)
both.fbok = true;
else
    both.fbok = false;
end

[cam1,cam2] = HSTcomputeParams(cam1,cam2,both);



cam1.figPos = [50 250 1310 560];%[50 50 1310 1100];
%[50 50 1310 560];%[50 50 1320 650];
cam1.figH = figure('pos',cam1.figPos,'PaperPositionMode','auto',...
    'visible','on');
cam2.figH = cam1.figH;



%% load 1D cut pixels Pix1 Pix2
% these come from HST/matlab/stereo/
if both.doFrameCuts
    addpath('precompute')
    load('HST1pix'),load('HST2pix')
    %preallocate
    cam1.keogram = zeros(cam1.nCol,both.nMutFrame,'uint16');
    cam2.keogram = zeros(cam2.nCol,both.nMutFrame,'uint16');
else
    Pix1 = []; Pix2 = []; %since we pass them into a fcn
end

%% setup HST1 panel
if both.playFrames

if cam1.use
%cam1.panelPos = [2 450 450 350];
% cam1.panelPos = [2 10 650 550];
%cam1.panelH = uipanel('parent',cam1.figH,'units','pixels','pos',cam1.panelPos);
cam1.panelH = uipanel('parent',cam1.figH,'units','normalized',...
                      'pos',[0.01 0.01 0.49 0.99]);

cam1.pbAx = axes('parent',cam1.panelH); %don't use nextplot yet, or image won't fill axis

cam1.imgH = imagesc(1:cam1.nCol,1:cam1.nRow,nan(cam1.nCol,cam1.nRow),...
            [cam1.minVal cam1.maxVal]);

%set(cam1.pbAx,'XTick',[],'ytick',[],... %erase tick marks
set(cam1.pbAx,'ydir','normal',...
              'nextplot','add') 

cam1.titleH = title('','interpreter','none');
cam1.cbH = colorbar('peer',cam1.pbAx);
set(get(cam1.cbH,'ylabel'),'string','HST1 16-bit data numbers [0..65535]')

%resize_fcn(cam1.panelH)

%plot 1D cut line
if both.doFrameCuts
%    plot(1:cam1.nCol,Pix1.LLSrow,'parent',cam1.pbAx) %6a
    plot(Pix1.xyAng(:,1),Pix1.xyAng(:,2),'parent',cam1.pbAx)
end

if both.plotPeaks
cam1.pPks = plot(nan,nan,'bo','markersize',10,'parent',cam1.pbAx); 

cam1.pTxt(1) = text(nan,nan,'1','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam1.pbAx);
cam1.pTxt(2) = text(nan,nan,'2','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam1.pbAx);
cam1.pTxt(3) = text(nan,nan,'3','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam1.pbAx);
end %if plotPeaks


end %if cam1.use

%% setup HST2 panel
if cam2.use
%     cam2.panelPos = [cam1.panelPos(1)+cam1.panelPos(3)+10 cam1.panelPos(2)...
%                      cam1.panelPos(3)                     cam1.panelPos(4)];
%cam2.panelH = uipanel('parent',cam1.figH,'units','pixels',...
%    'pos',cam2.panelPos);
cam2.panelH = uipanel('parent',cam1.figH,'units','normalized',...
                        'pos',[0.51 0.01 0.49 0.99]);

cam2.pbAx = axes('parent',cam2.panelH); 

cam2.imgH = imagesc(1:cam2.nCol,1:cam2.nRow,nan(cam2.nCol,cam2.nRow),...
            [cam2.minVal cam2.maxVal]);

%set(cam2.pbAx,'XTick',[],'ytick',[],... %erase tick marks
set(cam2.pbAx,'ydir','normal',...
              'nextplot','add')
cam2.titleH = title('','interpreter','none');
%axis('image')
cam2.cbH = colorbar('peer',cam2.pbAx);
set(get(cam2.cbH,'ylabel'),'string','HST2 14-bit data numbers [0..16384]')

%resize_fcn(cam2.panelH)


%plot 1D cut line
if both.doFrameCuts
   % plot(1:cam2.nCol,Pix2.LLSrow,'parent',cam2.pbAx)
    plot(Pix2.xyAng(:,1),Pix2.xyAng(:,2),'parent',cam2.pbAx)
end

if both.plotPeaks
cam2.pPks = plot(nan,nan,'bo','markersize',10,'parent',cam2.pbAx); 

cam2.pTxt(1) = text(nan,nan,'1','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam2.pbAx);
cam2.pTxt(2) = text(nan,nan,'2','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam2.pbAx);
cam2.pTxt(3) = text(nan,nan,'3','fontw','bold','HorizontalAlignment','center','VerticalAlignment','middle','parent',cam2.pbAx);
end %if plotPeaks

end %if cam2.use

%% setup manual transparancy layer

if both.manualPeaks
   cam1.hmPk = plot(nan,nan,'rs','markersize',5,'parent',cam1.pbAx);
   cam2.hmPk = plot(nan,nan,'rs','markersize',5,'parent',cam2.pbAx);
end


end %if playFrames
%% setup subtraction frame
if both.doFrameSub

both.substPanel = uipanel('parent',cam1.figH,'units','pixels',...
    'pos',[cam1.panelPos(1) cam1.panelPos(2)+(cam1.panelPos(4)+10)...
           cam1.panelPos(3) cam1.panelPos(4)]);
both.substAx = axes('parent',both.substPanel);

%camera frame sizes must be same!
both.substImgH = imagesc(1:cam1.nCol,1:cam1.nRow,nan(cam1.nCol,cam1.nRow),...
            [0 1]); %handles normalized values

both.substtH = title(both.substAx,'','interpreter','none');
%axis('image')
both.substcbH = colorbar(both.substAx);
set(get(both.substcbH,'ylabel'),'string','normalized intensity')

resize_fcn(both.substPanel)
end %if doFrameSub
%% setup Grad frame
if both.doFrameGrad
cam1 = makeGradPanels(cam1);
cam2 = makeGradPanels(cam2);
end %if do FrameGrad
%% colormap -----------------------------
colormap('gray') %set once, it affect whole figure

%% setup control panel 
% ctrl.panelPos = [2 610 800 120];
% ctrl.panelH = uipanel('parent',cam1.figH,'units','pixels','pos',ctrl.panelPos);
% ctrl.SkipPB = uicontrol('parent',ctrl.panelH,'style','pushbutton',...
%                 'Position',[5 5 50 18],'String','Skip','Callback',@skipBframe)

%% setup separate keogram figure
if both.doFrameCuts && ~both.keypressPlayback
    %overall keogram figure
   both.hfKeo = figure('Name','1D Keograms','pos',[20 20 1200 600]);
   
   %HST1 panel
   cam1.hpKeo = uipanel('parent',both.hfKeo,'units','pixels',...
       'pos',[1 1 600 500]);
   cam1.hKeoAx = axes('parent',cam1.hpKeo);
   cam1.hKeoImg = imagesc(1:cam1.nCol,uint32(cam1.firstFrameNum+cam1.pbInd-1),...
                          zeros(both.nMutFrame,cam1.nCol,'uint16'));
   set(cam1.hKeoAx,'ydir','normal')
   xlabel('1D cut pixel #')
   ylabel('frame #')
   title(['HST1 Keogram of 1D cut', datestr(both.tReqUT(1)),'UT to ',...
    datestr(both.tReqUT(end)),'UT'])

   cam1.hkeoCB = colorbar('peer',cam1.hKeoAx);
   set(get(cam1.hkeoCB,'ylabel'),'string','Data Numbers')
   
   %fix annoying exponential axis labels
   set(cam1.hKeoAx,'YTickLabel',sprintf('%7d|',get(cam1.hKeoAx,'YTick')))
   
   %HST2 panel
   cam2.hpKeo = uipanel('parent',both.hfKeo,'units','pixels',...
       'pos',[610 1 600 500]);
   cam2.hKeoAx = axes('parent',cam2.hpKeo);
   cam2.hKeoImg = imagesc(1:cam2.nCol,uint32(cam2.firstFrameNum+cam2.pbInd-1),...
                          zeros(both.nMutFrame,cam2.nCol,'uint16'));
   set(cam2.hKeoAx,'ydir','normal')
   title(['HST2 Keogram of 1D cut', datestr(both.tReqUT(1)),'UT to ',...
    datestr(both.tReqUT(end)),'UT'])
   xlabel('1D cut pixel #')
   ylabel('frame #')
   cam2.hkeoCB = colorbar('peer',cam2.hKeoAx);
   set(get(cam2.hkeoCB,'ylabel'),'string','Data Numbers')
   
   %fix annoying exponential axis labels
   set(cam2.hKeoAx,'YTickLabel',sprintf('%7d|',get(cam2.hKeoAx,'YTick')))
   
 %% setup stack plots
%{
 %overall keogram figure
   both.hfStack = figure('Name','1D Stack Plots','pos',[30 30 1200 600]);
   
   %HST1 panel
   cam1.hpStack = uipanel('parent',both.hfStack,'units','pixels',...
       'pos',[1 1 600 500]);
   cam1.hStackAx = axes('parent',cam1.hpStack);

   %HST2 panel
   cam2.hpStack = uipanel('parent',both.hfStack,'units','pixels',...
       'pos',[610 1 600 500]);
   cam2.hStackAx = axes('parent',cam2.hpStack);
%}
end




end


function cam = makeGradPanels(cam)
    cam.gradPanel = uipanel('parent',cam.figH,'units','pixels',...
    'pos',[cam.panelPos(1) cam.panelPos(2)+(cam.panelPos(4)+10)...
           cam.panelPos(3) cam.panelPos(4)]);
     cam.gradAx = axes('parent',cam.gradPanel);
     
cam.gradImgH = imagesc(1:cam.nCol,1:cam.nRow,nan(cam.nCol,cam.nRow),...
            [cam.minGrad cam.maxGrad]); %handles normalized values
%cam.gradImgH = subimage(1:cam.nCol,1:cam.nRow,nan(cam.nCol,cam.nRow));
%set(cam.gradImgH,'CDataMapping','scaled')

cam.gradtH = title(cam.gradAx,'Gradient','interpreter','none');

cam.gradcbH = colorbar('peer',cam.gradAx);
set(get(cam.gradcbH,'ylabel'),'string','normalized intensity')

end