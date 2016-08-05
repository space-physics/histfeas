.. image:: https://travis-ci.org/scienceopen/histfeas.svg?branch=master
    :target: https://travis-ci.org/scienceopen/histfeas

.. image:: https://codeclimate.com/github/scienceopen/histfeas/badges/gpa.svg
   :target: https://codeclimate.com/github/scienceopen/histfeas

.. image:: https://coveralls.io/repos/scienceopen/histfeas/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/scienceopen/histfeas?branch=master

================
hist-feasibility
================
Feasibility study for auroral tomography

:Author: Michael Hirsch
:Date: June 2015

This is a big program, contact me with any issues

.. contents::

.. image:: doc/montout.png
    :alt: montage of output

Installation
------------
::
   python setup.py develop

self-test::

    python tests/test.py

Usage notes
------------
if upon changing the ``in/*.ini`` files to make a new simulation, you get an error message
including
``use --ell command line option to save new Ell file``
rerun your simulation command, adding ``--ell`` to compute (one-time) the projection
matrix for your new simulation geometry.

Simulation Examples
-------------------

simulate flaming aurora with two cameras
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

  python RunHistfeas.py in/2cam_flame.ini out/test_flame2/ -m fwd optim png show h5

you can then look to the `Output Processing`_ section for how to load the HDF5 files
you just produced in ``out/test_flame2``

Real Data Examples
------------------

reading real data and displaying a live video
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

  python RunHistfeas.py in/apr14.ini out/apr14 -m realvid -a 0.1

reading real data and saving the joint image frames to disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

   python RunHistfeas.py in/apr14.ini out/apr14 -m realvid rawpng -a 0.1




Utility Examples
----------------

plot eigenprofiles from 2013 JGR and current transcar sim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

   python RunHistfeas.py in/jgr2013_2cam.ini /tmp -m eig eig1d -p  -f 0 1 1

   python RunHistfeas.py in/2cam_flame.ini /tmp -m eig eig1d -p --vlim 0 0 90 1000 1e-1 5e3 -f 0 1 1

Output selection (via -m command)
---------------------------------
combine the following commands as desired under the ``-m`` option to control the
type of program output

*Almost all of these* ``-m`` *options can be combined in various ways desired by the user*

**Simulation selection**

-m fwd      run foward model
-m optim    run optimization to estimate input quantities


**Graphics selection**

-m h5       dumps HDF5 files of the quantities selected in -m
-m eps      saves figures as eps
-m png      saves figures as png

**real data only**

-m realvid      both cameras in one figure
-m singleraw    each camera images individually, without axes (for powerpoint,posters, etc.)

**excitation rates plots**

-m eig         plot eigenprofiles
-m spectra     plot modeled auroral spectra modulated by the filter used.


Time selection
--------------

The simulation configuration in the in/\*.xlsx file may be very large. Maybe you want to pick
only a few times to run.

Example: to use only the first time step, use option ``-f 0 1 1`` which works like Python
``range()`` in selecting times from the spreadsheet Arc* tab.

Note on simulation time
~~~~~~~~~~~~~~~~~~~~~~~

The simulation time currently runs 10x faster than the columns in the in/\*.xlsx
under the Arc\* tabs. You should normally have the times of the Arc\* .xlsx columns
evenly spaced. If not, you can skip over the jump times by taking say every other time.

plot limit selection
--------------------
You may want to select fixed limits for your plots instead of the default autoscaling, particularly when
comparing a time series of plots.

====================================================  ===========
plot selection  description
====================================================  ===========
--vlim xmin xmax zmin zmax pmin pmax pmin1d pmax1d      limits for VER plots and eigenprofile plots (including 1-D)
--jlim min max min1d max1d                              x-axis limits for diff num flux plots (first two for 2-D, last two for 1-D)
--blim min max                                          intensity (y-axis) limits for brightness plots
====================================================  ===========

Plot explanation
-----------------
The plots you see under your out/ direction (assuming you used ``-m png`` or ``-m eps`` or the like)
follow this naming convention

``phifwd`` this is your "known" input differntial number flux of the electron precipitation
to the simulation (for real data, we don't have this)

``phiest`` this is the unobservable "unknown" we estimate with this program (for real and simulated data)

``pfwd`` and ``pest`` volume emission rate due to simulated / estimated flux respectively

``bfwd`` and ``best`` camera optical intensity due to simulated / estimated flux respectively

Our IEEE TGARS article (in review) details the math and algorithm.


Output Processing
------------------
HDF5 is the primary output along with PNGs of selected plots.
Some of the 1-D variables are duplicated because we don't know a-priori simulation parts will be run--disk space use is trivial, so we have left this alone.

The naming of the variables follows `Plot explanation`_

For Python, we have the hollow function ``loadAnalyze.py`` which loads the HDF5 data to call
the same ``analysehst.py`` that's used by the simulation online--good coding practice.

**Example of offline output processing**

::

  python loadAnalyze.py test/registration.h5



Calibration
-------------

1. ``rawDMCreader.py``  accesses the raw camera data and averages the selected frames and writes the average as a FITS file
2. The second line moves this FITS file to the user-selected calibration directory
3. The third line uses my wrapper and post-processing based on Astrometry.net to make an HDF5 file of the mapping from each pixel to sky coordinates (ra/dec and az/el).

**cam0**

::

   ./histutils/rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST0/2013-04-14T07-00-CamSer7196_frames_363000-1-369200.DMCdata -f 0 10 1 --avg --fits
   

   mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer7196_frames_363000-1-369200_mean_frames.fits ~/HST/calibration/hst0cal.fits
   
   ./astrometry/fits2azel.py -i ~/HST/calibration/hst0cal.fits --h5 -c 65.1186367 -147.432975 -t 2013-04-14T08:54:00Z --png


**cam1**

::

   ./histutils/rawDMCreader.py -i ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621.DMCdata -f 0 10 1 --avg --fits
   
   mv ~/HSTdata/DataField/2013-04-14/HST1/2013-04-14T07-00-CamSer1387_frames_205111-1-208621_mean_frames.fits ~/HST/calibration/hst1cal.fits
   
   ./astrometry/fits2azel.py -i ~/HST/calibration/hst1cal.fits --h5 -c 65.12657 -147.496908333 -t 2013-04-14T08:54:00Z --png


Install Errors
==============
Since the inherent function of this program is generating plots, you will need a windowing system::

	apt-get install libsm6 libxrender1 libfontconfig1

