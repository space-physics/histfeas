# hist-feasibility
Feasibility study for auroral tomography


 calibration:
-------------
cam0
```
 ./rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST0/2013-04-14T07-00-CamSer7196_frames_363000-1-369200.DMCdata -f 0 10 1 --avg --fits
 mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer7196_frames_363000-1-369200_mean_frames.fits ~/HST/calibration/hst0cal.fits
 /fits2azel.py -i ~/HST/calibration/hst0cal.fits --h5 -c 65.1186367 -147.432975 -t 2013-04-14T08:54:00Z --png
```

cam1
```
 ./rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621.DMCdata -f 0 10 1 --avg --fits
 mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621_mean_frames.fits ~/HST/calibration/hst1cal.fits
 ./fits2azel.py -i ~/HST/calibration/hst1cal.fits --h5 -c 65.12657 -147.496908333 -t 2013-04-14T08:54:00Z --png
```
