#!/usr/bin/env python
import h5py
from numpy import diff
from matplotlib.pyplot import figure,show,subplots
import seaborn as sns
sns.set_context('talk')

fn = '../out/085430-x-0-30/dump2013-04-14T08:54:30.226.h5'

cam = [{},{}]

with h5py.File(fn,'r') as f:
    cam[0]['angle'] = f['best']['angle'][0,:]
    cam[1]['angle'] = f['best']['angle'][1,:]


#%%
fg,axs = subplots(1,2)
ax = axs[0]
ax.set_title('camera angles')
for i,c in enumerate(cam):
    ax.plot(c['angle'],label='cam'+str(i))
ax.set_ylabel('angle [deg]')
ax.set_xlabel('pixel index #')
ax.autoscale(True,'x',True)
ax.legend()
#%%
ax = axs[1]
ax.set_title('diff(angle)')
for i,c in enumerate(cam):
    ax.plot(diff(c['angle']),label='cam'+str(i))
ax.set_ylabel('diff(angle) [deg]')
ax.set_xlabel('pixel index #')
ax.autoscale(True,'x',True)
ax.axhline(9/512,color='red',linewidth=2,linestyle='--',label='expected')
ax.legend()


show()