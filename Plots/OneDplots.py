#!/usr/bin/env python
import h5py
from numpy import diff,gradient,hypot
from matplotlib.pyplot import figure,show,subplots
import seaborn as sns
sns.set_context('talk')

def plotangles(fn):
#%% load data
    cam = [{},{}]

    with h5py.File(fn,'r') as f:
        cam[0]['angle'] = f['best']['angle'][0,:]
        cam[1]['angle'] = f['best']['angle'][1,:]
#%% plot angles
    fg,axs = subplots(1,2)
    ax = axs[0]
    ax.set_title('camera angles')
    for i,c in enumerate(cam):
        ax.plot(c['angle'],label='cam'+str(i))
    ax.set_ylabel('angle [deg]')
    ax.set_xlabel('pixel index #')
    ax.autoscale(True,'x',True)
    ax.legend()
#%% plot diff(angle)
    ax = axs[1]
    ax.set_title('diff(angle)')
    for i,c in enumerate(cam):
        ax.plot(diff(c['angle']),label='cam'+str(i))
    ax.set_ylabel('diff(angle) [deg]')
    ax.set_xlabel('pixel index #')
    ax.autoscale(True,'x',True)
    ax.axhline(9/512,color='red',linewidth=2,linestyle='--',label='expected')
    ax.legend()

def plotplatescale(flist):
#%% az, el
    fg,axs = subplots(len(flist),2,sharex=True,sharey=True)

    for i,c in enumerate(flist):
        fg.suptitle('azimuth & elevation')

        with h5py.File(c,'r') as f:
            fg.colorbar(axs[i,0].pcolormesh(f['az']),ax=axs[i,0]) #edgecolors too slow
            fg.colorbar(axs[i,1].pcolormesh(f['el']),ax=axs[i,1])
        axs[i,0].set_title('azimuth cam {}'.format(i))
        axs[i,1].set_title('elevation cam {}'.format(i))
    axs[0,0].autoscale(True,'both',True)
#%% diff(az, el)
    fg,axs = subplots(len(flist),2,sharex=True,sharey=True)

    for i,c in enumerate(flist):
        fg.suptitle('mag(gradient(azimuth & elevation))')

        with h5py.File(c,'r') as f:
            dAz = gradient(f['az'])
            dEl = gradient(f['el'])
        dazmag = hypot(dAz[0],dAz[1])
        delmag = hypot(dEl[0],dEl[1])

        fg.colorbar(axs[i,0].pcolormesh(dazmag,cmap='bwr'),ax=axs[i,0])
        fg.colorbar(axs[i,1].pcolormesh(delmag,cmap='bwr'),ax=axs[i,1])
        axs[i,0].set_title('azimuth cam {}'.format(i))
        axs[i,1].set_title('elevation cam {}'.format(i))
    axs[0,0].autoscale(True,'both',True)



if __name__ == '__main__':
    plotangles('../out/085430-x-0-30/dump2013-04-14T08:54:30.226.h5')
    plotplatescale(['../precompute/hst0cal.h5','../precompute/hst1cal.h5'])

    show()