#!/usr/bin/env python
from pathlib import Path
from wand.image import Image

def bfgsiterplot(plist,pfn,ofn):
    with Image() as anim:
        for p in plist:
            fn = p/pfn
            with Image(filename=str(fn)) as I:
                anim.sequence.append(I)

        print(f'saving {ofn}')         
        anim.save(filename=str(ofn))

def findplots(path,pat,pfn):
    path = Path(path).expanduser()
    plist = sorted(path.glob(pat+'*'))
    plist = [p for p in plist if p.is_dir()]# and (p/pfn).is_file()] # take only directories

    print(f'{len(plist)} files found in {path}')

    return plist

if __name__ == '__main__':
    path = '../out'
    pat = 'glitch0-'
    pfn = 'est2013-04-14T085430.226.png'
    
    plist = findplots(path,pat,pfn)
    print(plist)

    bfgsiterplot(plist,pfn,pat+'.gif')
