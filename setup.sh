#!/bin/sh

#may need this for tifffile to work
#sudo apt-get install libfreeimage3

conda install --file requirements.txt
pip install -r piprequirements.txt

# can be used in the future to refresh
git submodule foreach "(git checkout master; git pull)"

cd gridaurora/lowtran
make -f Makefile.f2py
