#!/bin/bash
set -u; set -e
# iterates over camera positions, using the step/impulse DNF input.

cam1=(0.5 1 1.5 2 2.5 3 3.5 4 5 6 7 8 9 10)
cam0=0

for c in ${cam1[*]}; do
    cx=$(printf "%0.2f" $c)
	./figure_impulse.py --cx $cam0 $cx -o ~/Dropbox/HirschThesis/impulse2/cam1_$cx
done


