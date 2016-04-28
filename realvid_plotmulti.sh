#!/bin/bash
set -e

# ./real_apr14.py -c realvid fwd png

#For a real video root directory, finds all the time results and makes multipage GIFs or TIFFs as you like from the PNGs.
# if you have very many PNGs, you might like to use FFV1 in an .avi instead of .gif.

rdir=$1

imstem=rawFrame*.png
fstem=fwd*.png

outfn=rawvid.gif


#declare -a flist=$(find $rdir -maxdepth 1 -mindepth 1 -name "$imstem"  | sort -V)
#declare -a glist=$(find $rdir -maxdepth 1 -mindepth 1 -name "$fstem"  | sort -V)
declare -a flist=$(ls $rdir/$imstem | sort -V)
declare -a glist=$(ls $rdir/$fstem  | sort -V)


for i in ${flist[*]}; do
echo $i
done

for i in {0..N-1}; do #assumes flist and glist are the same length
tmpfn=/tmp/$(printf "%03d" $i)tmp.png
echo "${flist[i]} ${glist[i]} -> $tmpfn"
convert -append "${flist[i]}" "${glist[i]}" tmpfn
done

convert -convert -delay 30 /tmp/*tmp.png $rdir/$outfn
rm /tmp/*tmp.png
