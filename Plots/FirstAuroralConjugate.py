#!/usr/bin/env python
"""
Willis 1996 http://adsabs.harvard.edu/abs/1996QJRAS..37..733W
"""
from numpy import arange,array
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import show,figure

sitella = [[-10.5, 122.8],
           [ 40.1, 117.4]]

Narclla = array([[41,118],
                 [41.5,123.5],
                 [40.5,127.5],
                 [38.5,132]])

Sarclla = array([[-26,  119],
                 [-26.5,123],
                 [-26,  126],
                 [-25,  129.5]
                 ])

ax= figure().gca()
m = Basemap(projection='merc',
              llcrnrlat=-40, urcrnrlat=45,
              llcrnrlon=90,urcrnrlon=160,
              lat_ts=20,
              resolution='l')
m.drawcoastlines()
m.drawcountries()
#m.drawmeridians(arange(0,360,30))
#m.drawparallels(arange(-90,90,30))
#%%
for lla in sitella:

    x,y = m(lla[1], lla[0])
    m.plot(x,y,'o',color='limegreen',markersize=12)
    m.plot
#%%
x,y = m(Narclla[:,1],Narclla[:,0])
m.plot(x,y,color='firebrick',linewidth=2)

x,y = m(Sarclla[:,1], Sarclla[:,0])
ax.plot(x,y,color='firebrick',linewidth=2)

ax.set_title('First Conjugate Auroral Observation 1770 CE')


show()