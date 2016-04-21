#!/bin/bash
set -e

#For a simulation root directory, finds all the camera position plot results and makes multipage GIFs or TIFFs as you like.

rdir=$1
n=$2

imstem=est$n.png


flist=$(find $rdir -maxdepth 2 -mindepth 2 -name $imstem | sort )
echo "converting ${#flist[@]} png under $rdir matching pattern $imstem to gif"
echo ${flist[*]}
#convert -delay 20 $flist $rdir/est$n.gif
