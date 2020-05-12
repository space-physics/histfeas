#!/usr/bin/env python
"""
Stormer's trigometric auroral height determination
Egeland 2013 p. 98

Explores effects of features at different heights vs. observer standoff distance.
Shows difficulty of fine-scale auroral observation:
    1. need close spaced caemras to unambiguously estimate fine scale (Semeter 2012)
    2. close cameras present difficult but solvable inversion problem (Hirsch 2015, 2017)
"""
import numpy as np
from matplotlib.pyplot import figure, show
import seaborn as sns

sns.set_context("talk", font_scale=1)
sns.set_style("whitegrid")
sns.set_palette("cubehelix")

b = 88  # angle to feature from observer #2
a = np.arange(45, b - 0.5, 0.01)  # degrees  # angle to feature from observer #1

d = np.array([3.4, 4.3, 20, 50])  # observer standoff distance [km]

br = np.radians(b)
ar = np.radians(a)
h = np.sin(ar) * np.sin(br) / np.sin(br - ar) * d[:, None]
h[(75 > h) | (h > 300)] = np.nan
#%%

fg = figure()
fg.clf()
axs = fg.subplots(2, 1, sharex=True)
ax = axs[0]
ax.plot(a, h.T)
ax.set_ylabel("altitude [km]")
ax.set_title("Størmer trigonometric altitude vs. observer standoff distance [km]")
ax.legend(d.astype(str), loc="best")

dh = np.diff(h)
ax = axs[1]
ax.plot(a[:-1], dh.T)
# ax.set_yscale('log')
# ax.set_ylim((0,100))
ax.set_ylabel(r"d$h$/d$\alpha$  [km/degree]")
ax.set_title("Angular sensitivity vs. observer standoff distance [km]")
ax.legend(d.astype(str), loc="best")
ax.set_xlabel(r"observer angle $\alpha$ [deg]")
show()

#%% symbolic
# Broken? Doesn't solve.
if False:
    from sympy import sin, solve, Symbol

    a = Symbol("a", real=True, positive=True)
    b = Symbol("b", real=True, positive=True)

    h = 100

    S = solve(sin(a) * sin(b) / sin(b - a) * d - h, (a, b))
    print(S)

#%% Human auditory reaction time
# I neglected to write down which article these came from
# was thinking about error due to human shutter reaction time to aurora/telephone stimulus
rt_ms = np.array([75, 85, 85, 85, 85, 85, 95, 95])
print(
    "Human auditory reaction time std() {:.1f} [ms]  mean: {:.1f}".format(
        rt_ms.std(), rt_ms.mean()
    )
)
