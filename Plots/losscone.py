#!/usr/bin/env python
"""
Loss cone angles vs. L-shell
"""
from matplotlib.pyplot import figure, show
from numpy import arange, arcsin, degrees
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('talk', font_scale=1.5)

L = arange(1.5, 8, 0.01)

alpha = arcsin((4*L**6 - 3*L**5)**(-1/4))


ax = figure().gca()
ax.plot(L, degrees(alpha), linestyle='--', color='black')
ax.fill_between(L, 0, degrees(alpha))
ax.set_title('loss cone angles vs. L-shell')
ax.set_xlabel('L-shell')
ax.set_ylabel('Loss cone [deg.]')
ax.set_xlim((L[0], L[-1]))
show()
