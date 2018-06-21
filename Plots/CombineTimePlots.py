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
from pathlib import Path
from wand.image import Image
from tempfile import NamedTemporaryFile
import subprocess
#from wand.exceptions import WandError
#from wand.display import display

MINHEIGHT = 500  # each element of canvas will be scaled up to at least this height--avoids widely disparate image sizes


def pngsmash(rdir, impat, fpat, outfn):
    rdir = Path(rdir).expanduser()
    outfn = Path(outfn).expanduser()

    if not outfn.suffix:
        raise ValueError('you must specify a suffix e.g.  .gif to the output filename so that Wand knows what format to write')

    ilist = sorted(rdir.glob(impat))
    flist = sorted(rdir.glob(fpat))

    assert len(ilist) == len(flist), 'unequal len() {} & {}    {} != {}'.format(ilist, flist, len(ilist), len(flist))

    with Image() as anim:
        for i, f in zip(ilist, flist):
            with Image(filename=str(i)) as I, Image(filename=str(f)) as F, Image() as J:
                # %% enforce minimum height (aspect-preserving scale increase if needed)
                for im in (I, F):
                    if im.height < MINHEIGHT:
                        im.transform(resize='x' + str(MINHEIGHT))
# %% compose composite canvas
                J.blank(max(I.width, F.width), I.height + F.height)  # blank canvas on which to place images
                J.composite(F, 0, 0)        # add data to canvas
                J.composite(I, 0, F.height)  # add video frame to canvas
                anim.sequence.append(J)

        if outfn.suffix == '.gif':
            print('writing', outfn)
            anim.save(filename=str(outfn))
        elif outfn.suffix in ('.avi', '.mp4', '.ogv', '.webm'):
            with NamedTemporaryFile(suffix='.gif') as f:  # forcing .gif temp since it's what Wand can handle
                print('using tempfile', f.name)
                anim.save(filename=f.name)  # the handle didn't work for some reason
                print('writing', outfn)
                # NOTE: -c:v ffv
                subprocess.call(['ffmpeg', '-i', f.name, '-c:v', 'ffv1', outfn])


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('rdir', help='directory where the sim/inversion output PNGs are')
    p.add_argument('ofn', help='output GIF filename')
    p.add_argument('-i', '--impat', help='glob pattern for raw image PNG', default='rawFrame*.png')
    p.add_argument('-f', '--fpat', help='glob pattern for inversion PNG', default='est*.png')
    p = p.parse_args()

    pngsmash(p.rdir, p.impat, p.fpat, p.ofn)


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
