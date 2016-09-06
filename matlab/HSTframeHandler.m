function [cam1,cam2,IntEl1,IntEl2] = HSTframeHandler(cam1,cam2,both,writerObj,t,Pix1,Pix2,IntEl1,IntEl2)

% t is loop iteration index
%% read next HST1/2 frames
[cam1Frame,cam1time,frame1Ind] = getFrameData(cam1,both,t,cam1.rotate,1);
[cam2Frame,cam2time,frame2Ind] = getFrameData(cam2,both,t,cam2.rotate,2);
% ======== put frame-by-frame code after this line
%% oblique xcor code goes here (Hanna)

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
[cam1,cam2,IntEl1,IntEl2]=intElTime(cam1,cam2,both,cam1Frame,cam2Frame,t,IntEl1,IntEl2);
%% (optional) plot peaks from 1D peak finder program
 plotPeaks(cam1,cam2,both,t)
%% let Matlab plots catch up
   pause(0.01) %arbitrary, used to let Matlab Java graphics "catch up"
%% video/tiff writing  
   videoHandler(cam1,cam2,both,writerObj,frame1Ind,frame2Ind)  
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
        %Pix_.LLSrow is the row to draw the pixel from for this column
        %--based on a linear least squares fit of the 1D line in the 
        % respective camera image planes intersecting
        % both sites and the magnetic zenith.
 %        cam1cut1D(i) = cam1Frame(Pix1.LLSrow(i),i);
 %        cam2cut1D(i) = cam2Frame(Pix2.LLSrow(i),i);
         cam1cut1D(i) = cam1Frame(Pix1.xyAng(i,2),Pix1.xyAng(i,1));
         cam2cut1D(i) = cam2Frame(Pix2.xyAng(i,2),Pix2.xyAng(i,1));
    end

    
%% make keogram of 1D cuts
% store each time instant in column vector so that x-dimension is "time"
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
    
%% user-specified image transforms    
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
%subtracts one frame from another (not really useful!!)
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