#!/usr/bin/env python
"""
iterates over camera positions, using the step/impulse differential number flux input.
"""
from histfeas import Path
import subprocess

rdir = 'out'
ini = '../in/2cam_split.ini'

rdir = Path(rdir).expanduser()
cam1 = (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10) # km
cam0 = 0 # km

for cx1 in cam1:
    cmd = ['../FigureMaker.py',ini,str(rdir/'cam1_{}'.format(cx1)),
                      '--cx',str(cam0),str(cx1), ]
    print(' '.join(cmd))
    subprocess.Popen(cmd)
