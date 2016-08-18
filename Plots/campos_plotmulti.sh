#!/bin/bash
set -u; set -e

#For a simulation root directory, finds all the camera position plot results and makes multipage GIFs or TIFFs as you like.

rdir=$1
n=$2

imstem=est$n.png


declare -a flist=$(find $rdir -maxdepth 2 -mindepth 2 -name "$imstem" | sort -V)
echo "converting png under $rdir matching pattern $imstem to $rdir/est$ngif"
#echo ${flist[*]}
convert -delay 30 $flist $rdir/est$n.gif
