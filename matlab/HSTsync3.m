function [cam1,cam2,both] = HSTsync3(cam1,cam2,both)


%cam.tUT is the UTC time each and every frame occurred

%day2sec = 86400; %maps seconds to datenum unit step

%setup UT times each frame occurred
if cam1.use
    if cam1.blessTimeShift %~=0
        warning(['HST1 Time Shifted by ',num2str(cam1.blessTimeShift),' seconds'])
    end
cam1.tUT = cam1.blessTimeShift/both.day2sec + ... %shifted to account for timing error--yikes!
            cam1.fullFileStart +...
            ((cam1.firstFrameNum:cam1.lastFrameNum)-1)*cam1.kineticSec/both.day2sec;

if ~isnan(cam1.tUT)
    display(['Camera 1 start/stop UTC: ',datestr(cam1.startUT),' ',datestr(cam1.stopUT)])
else %using synthetic file
end

else
   display('Camera 1 not used')
end
%% HST 2
if cam2.use
    if cam2.blessTimeShift %~=0
        warning(['HST2 Time Shifted by ',num2str(cam2.blessTimeShift),' seconds'])
    end
cam2.tUT =      cam2.blessTimeShift/both.day2sec + ... %shifted to account for timing error--yikes!        
                cam2.fullFileStart +...
                ((cam2.firstFrameNum:cam2.lastFrameNum)-1)*cam2.kineticSec/both.day2sec;

if ~isnan(cam2.tUT)
    display(['Camera 2 start/stop UTC: ',datestr(cam2.startUT),' ',datestr(cam2.stopUT)])
else %using synthetic file
end

else
    display('Camera 2 not used')
end
%determine mutual start/stop frame -- play only over UTC times that both
%sites have frames available
both.startUT = max([cam2.startUT cam1.startUT]); %who started last
both.stopUT = min([cam2.stopUT cam1.stopUT]); %who ended first

%make playback time steps -- based on the "simulated" UTC times that do not
%correspond exactly with either camera, necessarily.
both.tReqUT = both.startUT:both.simKineticSec/both.day2sec:both.stopUT;
both.nMutRawFrame = length(both.tReqUT);
both.fullFileIndicies = 1:cam1.nFrame;

if ~isnan(cam1.tUT)
display([int2str(both.nMutRawFrame),' Mutual frames available from ',datestr(both.startUT),'UT to ',...
    datestr(both.stopUT),'UT'])
else %using synthetic file
end

%adjust start/stop to user request
both.tReqUT(both.tReqUT<both.reqStartUT) = [];
both.tReqUT(both.tReqUT>both.reqStopUT) = [];
both.nReqFrame = length(both.tReqUT);

if ~isnan(cam1.tUT)
display(['Per user specification, displaying ',int2str(both.nReqFrame),'frames from '])
display([datestr(both.tReqUT(1)),'UT to ',...
    datestr(both.tReqUT(end)),'UT'])
else %synthetic, use first frame only 
end

%use nearest neighbor interpolation to find mutual frames to display--i.e.
%sometimes one camera will have frames repeated, while the other camera
%might skip some frames altogether

%cam1.tUT: UTC start exposure time of each frame
%cam1.pbInd: indicies of original data file corresponding to requested UTC
if cam1.use
    if ~isnan(cam1.tUT) %real data
        cam1.pbInd = interp1(cam1.tUT,1:cam1.nFrame,...
                               both.tReqUT,'nearest');
    else %synthetic
        cam1.pbInd = 1:cam1.nFrame;
    end
end

if cam2.use
    if ~isnan(cam1.tUT) %real data
        cam2.pbInd= interp1(cam2.tUT,1:cam2.nFrame,...
                            both.tReqUT,'nearest');
    else %synthetic
        cam2.pbInd = 1:cam2.nFrame;
    end
end



both.nMutFrame = length(cam1.pbInd); %both cam1 and cam2 will have same length after interp1

% display(['virtual playback rate user specified: ',num2str(both.pbFPS),'fps'])

display(['Playback kinetic time = ',num2str(both.simKineticSec),' sec.'])
display(['HST2 kinetic time = ',num2str(cam2.kineticSec),' sec.'])
display(['HST1 kinetic time = ',num2str(cam1.kineticSec),' sec.'])


end