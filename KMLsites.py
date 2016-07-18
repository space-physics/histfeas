"""
script to write instrument locations to KML for repeatable, extensible plots in Google Earth
"""
import simplekml
#
from histfeas import Path
from histfeas.io import campoint

kmlfn = 'cam.kml'
sites = {'HiST0':(65.1186367, -147.432975),
         'HiST1':(65.12657,   -147.496908333),
         'PFISR':(65.12992,   -147.47104)}

#%%
kmlfn = Path(kmlfn).expanduser()
kml = simplekml.Kml()

for name,latlon in sites.items():
    kml = campoint(kml,latlon,name,lkat=None)

print('saving {}'.format(kmlfn))
kml.save(str(kmlfn))