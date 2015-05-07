# hist-feasibility
Feasibility study for auroral tomography

This program should be runnable on any Python 2.7/3.4 platform. Primarily tested on Linux.

setup:
------------
```
pip install -r requirements.txt
git clone --recursive https://github.com/scienceopen/hist-feasbility.git
```

examples:
---------
draft 2015 SIMULATION commands were like:
```
python3 main_hist.py in/jgr2013_2cam_flame.xlsx /tmp --minev 150 -m fwd optim gfit eig eig1d ell eavg png --vlim -3.8 7.1 90 350 1e5 1e8 --jlim 0 7e4 --blim 0 2.5e9 -f 2 5 1 --ell
```

reading real data:
```
python3 main_hist.py in/jgr2013_realdata.xlsx /tmp -m realvid fwd optim png rawpng --vlim -3.8 7.1 90 350 0 30 --jlim nan 0.15 --blim 0 2500 -f 0 1 1
```

dump raw frames with time superimposed to disk without axes (for draft, -f 30 70 10)
```
python3 main_hist.py in/jgr2013_realdata.xlsx /tmp -m singleraw rawpng -f 30 70 10
```

plot eigenprofiles from 2013 JGR
```
python3 main_hist.py in/jgr2013_2cam.xlsx /tmp --minev 150 -m eig eig1d show -f 0 1 1
```

plotting hints:
---------------
```-m rawpng``` saves the real video frames you chose to PNG with annotations/axes

```-m realvid``` shows both cameras in one big figure

```-m singleraw``` shows each camera images individually, without axes (for powerpoint,posters, etc.)

 calibration:
-------------
The first program, rawDMCreader.py, accesses the raw camera data and averages the selected frames and writes the average as a FITS file

The second line moves this FITS file to the user-selected calibration directory

The third line uses my wrapper and post-processing based on Astrometry.net to make an HDF5 file of the mapping from each pixel to sky coordinates (ra/dec and az/el). 

cam0
```
 ./histutils/rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST0/2013-04-14T07-00-CamSer7196_frames_363000-1-369200.DMCdata -f 0 10 1 --avg --fits
 mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer7196_frames_363000-1-369200_mean_frames.fits ~/HST/calibration/hst0cal.fits
 ./astrometry/fits2azel.py -i ~/HST/calibration/hst0cal.fits --h5 -c 65.1186367 -147.432975 -t 2013-04-14T08:54:00Z --png
```

cam1
```
 ./histutils/rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621.DMCdata -f 0 10 1 --avg --fits
 mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621_mean_frames.fits ~/HST/calibration/hst1cal.fits
 ./astrometry/fits2azel.py -i ~/HST/calibration/hst1cal.fits --h5 -c 65.12657 -147.496908333 -t 2013-04-14T08:54:00Z --png
```

