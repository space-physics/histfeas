import simplekml
import logging
from matplotlib.pyplot import figure,clf
#
from pymap3d.coordconv3d import aer2geodetic#, aer2ecef

def planviewkml(cam,xKM,zKM,makeplot,figh,odir):
    """
    https://developers.google.com/kml/documentation/cameras
    """
    decimaterayfactor = 16

    az = 106.444916022 #TODO arbitrary should be from real camera pointing
    srange = 500e3 #[m]

    clr=['r','b','g','m']
    kclr = ['ff5c5ccd','ffff0000']
    ax = figure(figh).gca(); clf()

    kml1d = simplekml.Kml()
#%% setup camera (I preferred LookAt)
#        camview = skml.Camera(latitude=67.2,
#                              longitude=-147.2,
#                              altitude=100e3, #meters, int
#                              heading = 180., tilt = 10., roll = 0.,
#                              altitudemode = skml.AltitudeMode.relativetoground)
#%% setup LookAt
    lkat = simplekml.LookAt(latitude=65.111,
                           longitude=-147.465,
                           altitude=0,
                           heading=180,
                           range=4e3,
                           tilt=45)

    lla = []
    for C in cam:
      if C.usecam:
        #az is the temporary scalar defined above FIXME
        el = C.angle_deg[::decimaterayfactor] #double colon
        Np = el.size
        lla.append((C.lon,C.lat))
        # get LLA of pixel rays at 100km altitude
        latre,lonre,altre = aer2geodetic(az,el,srange,C.lat,C.lon,C.alt_m)
        # get ECEF of center ray at 90km (bottom of model space)
        #centazray = az #TODO
        #centelray = cam[ci].angle_deg[Np//2]
        #xrc,yrc,zrc = aer2ecef(centazray,centelray,zbottom,lat0,lon0,alt0)
#%% camera ground location
        kml1d = campoint(kml1d,(C.lat,C.lon),C.name,lkat)
#%% camera rays
        if 'kmlrays' in makeplot:
            for cri in range(Np):
                linestr = kml1d.newlinestring(name='')
                linestr.coords = [(C.lon,      C.lat,      C.alt_m),
                                  (lonre[cri], latre[cri], altre[cri])]
                linestr.altitudemode = simplekml.AltitudeMode.relativetoground
                linestr.style.linestyle.color = kclr[C.name]

        ax.plot(lonre,latre,'x',color=clr[C.name],markersize=6)
        ax.plot(C.lon,C.lat,'o',color=clr[C.name],markersize=12,label='cam{}'.format(C.name))

    ax.set_ylabel('WGS84 latitude [deg.]')
    ax.set_xlabel('WGS84 longitude [deg.]')
    ax.set_title('pixel locations at 100km altitude')
    ax.legend()
#%% setup line on ground connecting sites
    """
    https://developers.google.com/kml/faq#linestyle
    https://simplekml.readthedocs.org/en/latest/geometries.html#simplekml.LineString
    """
    ls = kml1d.newlinestring(name='3 km', coords=lla)
    ls.style.linestyle.width = 5
    ls.style.linestyle.color = simplekml.Color.yellow
    ls.style.labelstyle.scale= 2.5
    ls.style.labelstyle.color= simplekml.Color.white
    ls.style.labelstyle.gxlabelvisibility=1
    ls.visiblity=1
#%% write KML
    try:
        kmlfn = odir/'cam.kml'
        logging.info('saving {}'.format(kmlfn))
        kml1d.save(str(kmlfn))
    except Exception as e:
        logging.error('Error writing KML {}   {}'.format(kmlfn,e))

def campoint(kml,latlon,sitename='',lkat=None):
    """
    camera location points
    latlon: len=2 or 3 vector of WGS84 lat,lon
    """
    bpnt = kml.newpoint(name=sitename,
                          #description = 'camera {} location'.format(C.name),
                        coords = [(latlon[1],latlon[0])],
                        altitudemode = simplekml.AltitudeMode.clamptoground)
    bpnt.style.iconstyle.icon.href='http://maps.google.com/mapfiles/kml/shapes/arrow.png'
    # 'http://maps.google.com/mapfiles/kml/paddle/pink-blank.png'
    bpnt.style.iconstyle.scale = 2.0
    bpnt.style.labelstyle.size= 2.5

    #bpnt.camera = camview
    if lkat is not None:
        bpnt.lookat = lkat

    return kml