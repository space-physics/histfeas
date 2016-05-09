#!/bin/env python3
"""
./real_apr14.py -c realvid fwd png

For a real video root directory, finds all the time results and makes multipage GIFs or TIFFs as you like from the PNGs.
if you have very many PNGs, you might like to use FFV1 in an .avi instead of .gif.
"""
from pathlib import Path
from sys import argv
from wand.image import Image
from wand.display import display

rdir=Path(argv[1])

imglob='rawFrame*.png'
fglob='fwd*.png'

outfn='rawvid.gif'

ilist = list(rdir.glob(imglob))
flist = list(rdir.glob(fglob))

assert len(ilist)==len(flist), 'unequal {} and fglob {} lengths: {}  {}'.format(ilist,flist,len(ilist),len(flist))

ilist.sort(); flist.sort()

with Image() as anim:
    for i,f in zip(ilist,flist):
       with Image(filename=str(i)) as I, Image(filename=str(f)) as F, Image() as J:
           J.blank(max(I.width,F.width), I.height + F.height)
           J.composite(F,0,0)
           J.composite(I,0,F.height)
           anim.sequence.append(J)

    print('writing {}'.format(outfn))
    anim.save(filename=outfn)



#tmpfn=/tmp/$(printf "%03d" $i)tmp.png
#echo "${flist[i]} ${glist[i]} -> $tmpfn"
#convert -append "${flist[i]}" "${glist[i]}" tmpfn
#done

#convert -convert -delay 30 /tmp/*tmp.png $rdir/$outfn
#rm /tmp/*tmp.png
