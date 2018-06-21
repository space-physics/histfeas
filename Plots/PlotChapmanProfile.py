#!/usr/bin/env python
"""
simple plot of Chapman profile
Michael Hirsch
"""
from numpy import arange
from matplotlib.pyplot import subplots, show
import seaborn as sns
sns.set_context('talk')
#
from histfeas.arcclass import ChapmanArc


def plotchapman(Wkm, H, X0, Z0):
    x = arange(-5, 5, 0.05)
    z = arange(60, 300, 5)

    Ne = ChapmanArc(Wkm, H, X0, Z0, xKM=x, zKM=z, xshape='gaussian')

    fg, axs = subplots(1, 2, sharey=True)

    ax = axs[0]
    h = ax.pcolormesh(x, z, Ne)
    fg.suptitle('example Chapman profile', fontsize='x-large')
    ax.set_title('2-D vertical slice')
    ax.set_xlabel('x [km]')
    ax.set_ylabel('z [km]')
    c = fg.colorbar(h, ax=ax)
    c.set_label('Normalized')
    ax.autoscale(True, tight=True)

    ax = axs[1]
    ax.plot(Ne[:, 99], z)
    ax.set_title('1-D vertical cut')
    #ax.set_ylabel('z [km]')
    ax.set_xlabel('normalized')


if __name__ == '__main__':
    plotchapman(Wkm=0.5, H=20, X0=0, Z0=140)
    show()
