#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots first recorded auroral conjugate plausible geometry.

Willis 1996 http://adsabs.harvard.edu/abs/1996QJRAS..37..733W

FIXME: add geomagnetic coordinates for 1770

Chinese font: apt install fonts-wqy-zenhei
"""
import numpy as np
import cartopy

#
# import matplotlib as mpl
# font_name = "WenQuanYi Zen Hei"
# mpl.rcParams['font.family']=font_name
#
# import matplotlib.font_manager as mfm
# ch_font = mfm.FontProperties(fname="/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc")
#
# from matplotlib import font_manager
# fontP = font_manager.FontProperties()
# fontP.set_family('WenQuanYi Zen Hei')
#
from matplotlib.pyplot import show, figure

#
from pymap3d import geodetic2aer

#
PROJ = cartopy.crs.PlateCarree()  # arbitrary


def main():
    # %% from paper
    sitella = {"HMS Endeavour": [-10.45, 122.82], "冀县 Ji-zhou": [40.1, 117.4]}

    Narclla = np.array([[41, 118], [41.5, 123.5], [40.5, 127.5], [38.5, 132]])

    Sarclla = np.array([[-26, 119], [-26.5, 123], [-26, 126], [-25, 129.5]])
    # %% simple calculations for presumable LOS
    archeight = 350e3  # m, assumed
    aer = geodetic2aer(Sarclla[1, 0], Sarclla[1, 1], archeight, *sitella["HMS Endeavour"], 0.0)
    print(f"aurora elevation angle (deg)  {aer[1]:.1f}")
    # %%
    ax = figure().gca(projection=PROJ)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=":")

    ax.set_extent((90, 160, -40, 45))

    #%% sites
    for o, lla in sitella.items():
        ax.plot(*lla[::-1], "o", color="limegreen", markersize=12, transform=PROJ)
        ax.annotate(
            o,
            xy=lla[::-1],
            xytext=(3, 3),
            textcoords="offset points",
            ha="right",
            family="WenQuanYi Wen Hei",
        )
    #%% aurora
    ax.plot(Narclla[:, 1], Narclla[:, 0], color="firebrick", linewidth=2, transform=PROJ)

    ax.plot(Sarclla[:, 1], Sarclla[:, 0], color="firebrick", linewidth=2, transform=PROJ)

    ax.set_title("First Conjugate Auroral Observation 1770 CE")
    # fontproperties=fontP)


if __name__ == "__main__":

    main()

    show()
