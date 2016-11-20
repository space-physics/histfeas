#!/usr/bin/env python
import subprocess
from setuptools import setup
from os import chdir

try:
    subprocess.call(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    pass
#%%
# FIXME: trick is to have them in order from no prereq to full prereq
for p in ['https://github.com/scienceopen/themisasi',
          'https://github.com/scienceopen/histutils',
        'https://github.com/scienceopen/lowtran',
        'https://github.com/scienceopen/pymap3d',
        'https://github.com/scienceopen/astrometry_azel',
        'https://github.com/scienceopen/CVutils',
        'https://github.com/scienceopen/gridaurora',
        'https://github.com/scienceopen/transcarread',
        'https://github.com/scienceopen/gaussfitter',
        'https://github.com/scienceopen/pyimagevideo',
        'https://github.com/scienceopen/dascutils',
        'https://github.com/scienceopen/pybashutils']:

    chdir('..')
    p.rstrip('/')
    subprocess.call(['git','clone', p])

    chdir(p.split('/')[-1])
    subprocess.call(['git','pull']) #in case it was already installed

    subprocess.call(['python','setup.py','develop']) #FIXME is there an API to do this?

chdir('../histfeas')
#%%
setup(name='histfeas',
      packages=['histfeas'],
	  description='Feasibility study for HiST auroral tomography system',
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/histfeas',
	  install_requires=['Wand','pathvalidate','geopy','simplekml'],
	  extras_require = {'tifffile':['tifffile']},
      dependency_links = [ ],
	  )

