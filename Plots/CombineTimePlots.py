#!/usr/bin/env python
"""
example creation of PNGs for use with this program:
./FigureMaker.py in/apr14T085454.ini

For a real video root directory, finds all the time results and makes multipage GIFs or TIFFs as you like from the PNGs
if you have very many PNGs, you might like to use FFV1 in an .avi instead of .gif.

Motivation: I have two separate programs (histfeas, histutils) that each create complicated figures. There is NOT currently
http://stackoverflow.com/questions/22521560/how-to-combine-several-matplotlib-figures-into-one-figure
a good way to combine two figures into one, even by grabbing axes.

Recommendation was to do this--create PNGs and smash together in post-processing.
"""
from histfeas import Path
from wand.image import Image
#from wand.display import display

def pngsmash(rdir,impat,fpat,outfn):
    rdir = Path(rdir).expanduser()

    ilist = sorted(rdir.glob(impat))
    flist = sorted(rdir.glob(fpat))

    assert len(ilist)==len(flist), 'unequal len() {} & {}    {} != {}'.format(ilist,flist,len(ilist),len(flist))

    with Image() as anim:
        for i,f in zip(ilist,flist):
           with Image(filename=str(i)) as I, Image(filename=str(f)) as F, Image() as J:
               J.blank(max(I.width,F.width), I.height + F.height)
               J.composite(F,0,0)
               J.composite(I,0,F.height)
               anim.sequence.append(J)

        print('writing {}'.format(outfn))
        anim.save(filename=outfn)

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('rdir',help='directory where the sim/inversion output PNGs are')
    p.add_argument('ofn',help='output GIF filename')
    p.add_argument('-i','--impat',help='glob pattern for raw image PNG',default='rawFrame*.png')
    p.add_argument('-f','--fpat',help='glob pattern for inversion PNG',default='est*.png')
    p = p.parse_args()

    pngsmash(p.rdir,p.impat,p.fpat,p.ofn)


"""
old way
#!/bin/sh

tmpfn=/tmp/$(printf "%03d" $i)tmp.png
echo "${flist[i]} ${glist[i]} -> $tmpfn"
convert -append "${flist[i]}" "${glist[i]}" tmpfn
done

convert -convert -delay 30 /tmp/*tmp.png $rdir/$outfn
rm /tmp/*tmp.png
"""