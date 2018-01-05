#!/usr/bin/env python
"""
Plots first recorded auroral conjugate plausible geometry.

Willis 1996 http://adsabs.harvard.edu/abs/1996QJRAS..37..733W
"""
import numpy as np
import cartopy
from matplotlib.pyplot import show,figure
#
from pymap3d import geodetic2aer
#
PROJ = cartopy.crs.PlateCarree()  # arbitrary

def main():
# %% from paper
    sitella = np.array([[-10.5, 122.8],
                        [ 40.1, 117.4]])

    Narclla = np.array([[41,118],
                        [41.5,123.5],
                        [40.5,127.5],
                        [38.5,132]])

    Sarclla = np.array([[-26,  119],
                        [-26.5,123],
                        [-26,  126],
                        [-25,  129.5]
                     ])
# %% simple calculations for presumable LOS
    archeight = 300e3 # m, assumed
    aer = geodetic2aer(Sarclla[1,0], Sarclla[1,1], archeight,
                       sitella[0,0], sitella[0,1], 0.)
    print(f'aurora elevation angle (deg)  {aer[1]:.1f}')
# %%
    ax= figure().gca(projection=PROJ)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')

    ax.set_extent((90,160,-40,45))

    #%% sites
    for lla in sitella:
        ax.plot(lla[1], lla[0],'o',
                color='limegreen',markersize=12,
                transform=PROJ)
    #%% aurora
    ax.plot(Narclla[:,1], Narclla[:,0],
           color='firebrick',linewidth=2,
           transform=PROJ)

    ax.plot(Sarclla[:,1], Sarclla[:,0],
            color='firebrick',linewidth=2,
            transform=PROJ)

    ax.set_title('First Conjugate Auroral Observation 1770 CE')

if __name__ == '__main__':

    main()

    show()