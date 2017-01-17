#!/usr/bin/env python
"""
Stormer's trigometric auroral height determination
Egeland 2013 p. 98
"""
from sympy import sin,solve,Symbol
import numpy as np
from matplotlib.pyplot import figure,show

a = Symbol("a",real=True,positive=True)
#b = Symbol("b",real=True,positive=True)
b = np.pi/2-0.1

d = 25
h=100

S=solve(sin(a) * sin(b) / sin(b-a) * d - h,(a,b) )
print(S)