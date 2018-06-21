#!/usr/bin/env python
import h5py
from numpy import diff, gradient, hypot, meshgrid, arange
from matplotlib.pyplot import figure, show, subplots
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set_context('talk')


def plotangles(fn):
    # %% load data
    cam = [{}, {}]

    with h5py.File(fn, 'r') as f:
        cam[0]['angle'] = f['best']['angle'][0, :]
        cam[1]['angle'] = f['best']['angle'][1, :]
# %% plot angles
    fg, axs = subplots(1, 2)
    ax = axs[0]
    ax.set_title('camera angles')
    for i, c in enumerate(cam):
        ax.plot(c['angle'], label='cam'+str(i))
    ax.set_ylabel('angle [deg]')
    ax.set_xlabel('pixel index #')
    ax.autoscale(True, 'x', True)
    ax.legend()
# %% plot diff(angle)
    ax = axs[1]
    ax.set_title('diff(angle)')
    for i, c in enumerate(cam):
        ax.plot(diff(c['angle']), label='cam'+str(i))
    ax.set_ylabel('diff(angle) [deg]')
    ax.set_xlabel('pixel index #')
    ax.autoscale(True, 'x', True)
    ax.axhline(9/512, color='red', linewidth=2, linestyle='--', label='expected')
    ax.legend()


def plotplatescale(flist, ptype='mesh'):
    # %% az, el
    if False:
        fg, axs = subplots(len(flist), 2, sharex=True, sharey=True)

        for i, c in enumerate(flist):
            fg.suptitle('azimuth & elevation')

            with h5py.File(c, 'r') as f:
                fg.colorbar(axs[i, 0].pcolormesh(f['az']), ax=axs[i, 0])  # edgecolors too slow
                fg.colorbar(axs[i, 1].pcolormesh(f['el']), ax=axs[i, 1])
            axs[i, 0].set_title('azimuth cam {}'.format(i))
            axs[i, 1].set_title('elevation cam {}'.format(i))
        axs[0, 0].autoscale(True, 'both', True)
# %% diff(az, el)
    azmesh(flist, ('az', 'el'))
# %% diff(ra,dec)
    azmesh(flist, ('ra', 'dec'))


def azmesh(flist, names, ptype='mesh'):

    if ptype == 'pcolor':
        fg, axs = subplots(len(flist), 2, sharex=True, sharey=True, projection='3d')
    elif ptype == 'mesh':
        fg = figure()
    j = 1  # for mesh
    for i, c in enumerate(flist):
        fg.suptitle('mag(gradient({}))'.format(names))

        with h5py.File(c, 'r') as f:
            dU = gradient(f[names[0]])
            dV = gradient(f[names[1]])
        dUmag = hypot(dU[0], dU[1])
        dVmag = hypot(dV[0], dV[1])

        if ptype == 'pcolor':
            fg.colorbar(axs[i, 0].pcolormesh(dUmag, cmap='bwr'), ax=axs[i, 0])
            fg.colorbar(axs[i, 1].pcolormesh(dVmag, cmap='bwr'), ax=axs[i, 1])
            axs[i, 0].set_title('RA cam {}'.format(i))
            axs[i, 1].set_title('DEC cam {}'.format(i))
            axs[i, 0].autoscale(True, 'both', True)
        elif ptype == 'mesh':
            X, Y = meshgrid(arange(dUmag.shape[1]), arange(dUmag.shape[0]))

            ax = fg.add_subplot(len(flist), 2, j, projection='3d')
            ax.plot_surface(X, Y, dUmag)
            ax.set_title('{} cam {}'.format(names[0], i))
            j += 1

            ax = fg.add_subplot(len(flist), 2, j, projection='3d')
            ax.plot_surface(X, Y, dVmag)
            ax.set_title('{} cam {}'.format(names[1], i))
            j += 1


if __name__ == '__main__':
    plotangles('../out/085430-x-0-30/dump2013-04-14T08:54:30.226.h5')
    plotplatescale(['../precompute/hst0cal.h5', '../precompute/hst1cal.h5'])

    show()
