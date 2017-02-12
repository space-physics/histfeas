#!/usr/bin/env python
import subprocess # need for git
from pathlib import Path
from setuptools import setup

req = ['Wand','pathvalidate','geopy','simplekml',
       'pymap3d','sciencedates','histutils','astrometry_azel','morecvutils','gridaurora',
       'nose','numpy','h5py','scipy','pandas','xarray','matplotlib','seaborn','astropy']
# leave astropy in here for gaussfitter


#%%
# FIXME: trick is to have them in order from no prereq to full prereq
for p in [
        'https://github.com/scienceopen/pybashutils',
        'https://github.com/scienceopen/pyimagevideo',
          'https://github.com/scienceopen/lowtran',          'https://github.com/scienceopen/dmcutils',
          'https://github.com/scienceopen/themisasi',
        'https://github.com/scienceopen/transcarread',
        'https://github.com/scienceopen/gaussfitter',
        'https://github.com/scienceopen/dascutils',]:

    cwd = Path('..') / p.split('/')[-1]
    print('\n {} \n'.format(cwd))

    if not cwd.is_dir():
        subprocess.run(['git','clone', p],cwd='..')

    subprocess.call(['git','pull'], cwd=str(cwd)) #in case it was already installed

    subprocess.call(['python','setup.py','develop'], cwd=str(cwd)) #FIXME is there an API to do this?

#%%
setup(name='histfeas',
      packages=['histfeas'],
	  description='Feasibility study for HiST auroral tomography system',
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/histfeas',
	  install_requires=req,
	  extras_require = {'tifffile':['tifffile']},
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',
      ],
	  )


