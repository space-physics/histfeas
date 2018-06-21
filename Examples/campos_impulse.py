#!/usr/bin/env python
"""
iterates over camera positions, using the step/impulse differential number flux input.
"""
import multiprocessing
import concurrent.futures
import subprocess
from pathlib import Path

rdir = 'out'
ini = '../in/2cam_split.ini'
cam1 = (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10)  # km
cam0 = 0  # km
# %%
Ncpu = multiprocessing.cpu_count()//2  # //2 makes one thread per CPU for 2 thread Intel Hyperthreading
rdir = Path(rdir).expanduser()

print(f'using {Ncpu} CPU')


def runhist(cx1):
    cmd = ['../FigureMaker.py', ini, str(rdir/f'cam1_{cx1}'),
           '--cx', str(cam0), str(cx1)]
    print(' '.join(cmd))
    subprocess.run(cmd)


def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=Ncpu) as executor:
        executor.map(runhist, cam1)


if __name__ == '__main__':
    main()
