#!/usr/bin/env python
import subprocess # need for git
from pathlib import Path
from setuptools import setup

try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)
    import pip
    pip.main(['install','-r','requirements.txt'])


#%%
# FIXME: trick is to have them in order from no prereq to full prereq
for p in [
        'https://github.com/scienceopen/pybashutils'
        'https://github.com/scienceopen/pyimagevideo',
          'https://github.com/scienceopen/lowtran',
        'https://github.com/scienceopen/CVutils',
        'https://github.com/scienceopen/astrometry_azel',
          'https://github.com/scienceopen/histutils',
          'https://github.com/scienceopen/dmcutils',
          'https://github.com/scienceopen/themisasi',
        'https://github.com/scienceopen/gridaurora',
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
	  install_requires=['Wand','pathvalidate','geopy','simplekml',
                      'pymap3d'],
	  extras_require = {'tifffile':['tifffile']},
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',
      ],
	  )


