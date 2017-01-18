#!/usr/bin/env python
"""
Stormer's trigometric auroral height determination
Egeland 2013 p. 98
"""
from sympy import sin,solve,Symbol
import numpy as np
from matplotlib.pyplot import figure,subplots,show
import seaborn as sns
sns.set_context('talk',font_scale=1.5)
sns.set_style('whitegrid')
sns.set_palette('cubehelix')

b = 88
a = np.arange(45,b-.5,0.01) #degress


d = np.array([3.4,4.3,20,50])

b = np.radians(b); a = np.radians(a)
h = np.sin(a) * np.sin(b) / np.sin(b-a) * d[:,None]
h[(75>h) | (h>300)] = np.nan
#%%

fg=figure(1); fg.clf()
ax = fg.gca()
ax.plot(np.degrees(a),h.T)
#ax.set_ylim((75,500))
ax.set_ylabel('altitude [km]')
ax.set_xlabel(r'observer angle $\alpha$ [deg]')


fg=figure(2); fg.clf()
fg,axs = subplots(2,1,num=2,sharex=True)
ax=axs[0]
ax.plot(np.degrees(a),h.T)
ax.set_ylim((70,310))
ax.set_xlim((60,None))
ax.set_ylabel('altitude [km]')
ax.set_title('Størmer trigonometric altitude estimation')
ax.legend(d.astype(str),loc='best')

dh = np.diff(h)
ax=axs[1]
ax.plot(np.degrees(a[:-1]),dh.T)
#ax.set_yscale('log')
#ax.set_ylim((0,100))
ax.set_ylabel(r'd$h$/d$\alpha$')
ax.set_title('Angular sensitivity vs. separation distance [km]')
ax.legend(d.astype(str),loc='best')
ax.set_xlabel(r'observer angle $\alpha$ [deg]')
show()

#%% symbolic
if 0:
    a = Symbol("a",real=True,positive=True)
    b = Symbol("b",real=True,positive=True)

    h=100
    
    S=solve(sin(a) * sin(b) / sin(b-a) * d - h,(a,b) )
    print(S)
    
#%% reaction time
rt_ms = np.array([75,85,85,85,85,85,95,95])
print('reaction time std() {:.1f} [ms]  mean: {:.1f}'.format(rt_ms.std(),rt_ms.mean()))
