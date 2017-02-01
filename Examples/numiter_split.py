#!/usr/bin/env python
"""
iterates over number of data inversion iterations, using the single splitting arc differential number flux input.
"""
import multiprocessing
import concurrent.futures
import subprocess
from pathlib import Path

load=False

rdir = 'out'
ini = '../in/2cam_split.ini'
outpat = 'bfgs_iter{}'

Niter = list(range(10)) + list(range(10,20,5)) + list(range(20,50,10)) + list(range(50,100,20))
#%%
Ncpu=multiprocessing.cpu_count()//2  # //2 makes one thread per CPU for 2 thread Intel Hyperthreading
rdir = Path(rdir).expanduser()
ini = Path(ini).expanduser()
assert ini.is_file()

print('using {} CPU'.format(Ncpu))

def runhist(nit):
    """
    This takes minutes to compute, and results are saved to HDF5.
    Note errors may be silent--try running without "concurrent" i.e. uncomment line 49
    """
    cmd = ['../FigureMaker.py', str(ini), str(rdir / outpat.format(nit)),'-m','fwd','optim',
           '--iter',str(nit)]
    print(' '.join(cmd))
    subprocess.run(cmd)

def loadhist(nit):
    """
    long-computation time results are loaded from HDF5, for faster post-processing
    """
    cmd = ['../FigureMaker.py', str((rdir / outpat.format(nit))/ini.name),'-m','fwd','optim',
           '--load']

    print(' '.join(cmd))
    subprocess.run(cmd)

def main(load):
#%% compute
    if not load:
        print("RUNNING COMPUTATIONS")
        #runhist(Niter[0])
        with concurrent.futures.ProcessPoolExecutor(max_workers=Ncpu) as executor:
            executor.map(runhist, Niter)
    else:
        print('SKIPPING COMPUTATIONS')
#%% analyze
    #odirs = [str(rdir / outpat.format(n)) for n in Niter]
    for nit in Niter:
        loadhist(nit)

if __name__ == '__main__':
    main(load)
