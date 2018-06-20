#!/usr/bin/env python
"""
Energy estimation from ASI such as DASC and THEMIS
using xarray as base class instead of custom classes.

./EnergyASI.py in/ASI-pkr-gako.ini
"""
import matplotlib
matplotlib.use('Agg')
print(matplotlib.get_backend())
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette("cubehelix")
sns.set(context='paper', style='whitegrid',font_scale=2,
        rc={'image.cmap': 'cubehelix_r'})
#
from histfeas import userinput
from histfeas.main_hist import doSim


if __name__ == '__main__':
    P = userinput()

    doSim(P)

    show()
