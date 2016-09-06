function [HST1,HST2,Targ,TargAlt] = LSEtarget2(HST1,HST2,Targ,both,t)
r2d = 180/pi;

%make hypothetic endpoint for auroral lines of view

[HST1.tx,HST1.ty,HST1.tz] = aer2ecef(Targ.Az1,Targ.El1,HST1.MaxAuroralSlantRange,...
                                     HST1.lla(1),HST1.lla(2),HST1.lla(3),...
                                     referenceEllipsoid('wgs84'),'degrees');
                                 
[HST2.tx,HST2.ty,HST2.tz] = aer2ecef(Targ.Az2,Targ.El2,HST2.MaxAuroralSlantRange,...
                                     HST2.lla(1),HST2.lla(2),HST2.lla(3),...
                                     referenceEllipsoid('wgs84'),'degrees');
                                 
if both.diagLSE
   [HST1.tLLA(1),HST1.tLLA(2),HST1.tLLA(3)] = ecef2geodetic(HST1.tx,HST1.ty,HST1.tz, referenceEllipsoid('wgs84'));
   HST1.tLLA(1:2) = HST1.tLLA(1:2)*r2d;
   
   [HST2.tLLA(1),HST2.tLLA(2),HST2.tLLA(3)] = ecef2geodetic(HST2.tx,HST2.ty,HST2.tz, referenceEllipsoid('wgs84'));
   HST2.tLLA(1:2) = HST2.tLLA(1:2)*r2d;
   
kmlwriteline([HST1.fn(2:11),'t',int2str(t),'ray1.kml'],[HST1.tLLA(1),HST1.lla(1)],...
                                                       [HST1.tLLA(2),HST1.lla(2)],...
                                                       [HST1.tLLA(3),HST1.lla(3)],...
     'AltitudeMode','RelativeToGround','Name',['HST1t',int2str(t)])
kmlwriteline([HST2.fn(2:11),'t',int2str(t),'ray2.kml'],[HST2.tLLA(1),HST2.lla(1)],...
                                                       [HST2.tLLA(2),HST2.lla(2)],...
                                                       [HST2.tLLA(3),HST2.lla(3)],...
     'AltitudeMode','RelativeToGround','Name',['HST2t',int2str(t)])
end
                        
% for all the pairwise HST1,HST2 observations, estimate the ECEF of the
% target
NpixIter = length(HST1.tx); AlotofPix = NpixIter>10;
for i = 1:NpixIter

[Pint]= lineIntersect3D([HST1.x,HST1.y,HST1.z;...
                              HST2.x,HST2.y,HST2.z],...
                             [HST1.tx(i),HST1.ty(i),HST1.tz(i);...
                              HST2.tx(i),HST2.ty(i),HST2.tz(i)]);
%{  
    [Xp,Yp,dp] = line2closest([HST1.x,HST1.y,HST1.z],[HST1.tx(i),HST1.ty(i),HST1.tz(i)],...
                        [HST2.x,HST2.y,HST2.z],[HST2.tx(i),HST2.ty(i),HST2.tz(i)]);
%}
[Targ.lla(i,1),Targ.lla(i,2),Targ.lla(i,3)] = ecef2geodetic(Pint(1),Pint(2),Pint(3),...
                                        referenceEllipsoid('wgs84'));
Targ.lla(i,1:2) = Targ.lla(i,1:2).*r2d;
Targ.x(i)=Pint(1); Targ.y(i)=Pint(2); Targ.z(i)=Pint(3);

Hra.kmlName{i} = int2str(i);%[num2str(xPix1(i)),',',num2str(yPix1(i))];

if ~AlotofPix    
display('------------------------------------------------------------')
display(['Target #',int2str(i),' Targ.Az1=',num2str(Targ.Az1(i)),' deg. and Targ.El1=',num2str(Targ.El1(i)),' deg.'])
display(['and Targ.Az2=',num2str(Targ.Az2(i)),' deg. and Targ.El2=',num2str(Targ.El2(i)),' deg.'])
display(['WGS84 Lat,Lon of auroral feature is [',num2str(Targ.lla(i,1),'%0.2f'),',',...
    num2str(Targ.lla(i,2),'%0.2f'),'], WGS84 alt=',num2str(Targ.lla(i,3)/1e3,'%0.2f'),' km.'])
end

end
TargAlt = Targ.lla(:,3);
end