#!/usr/bin/env python
install_requires=['numpy','h5py','scipy','pandas','xarray','matplotlib', 'seaborn', 'astropy',
     'Wand','pathvalidate','geopy','simplekml',
     'pymap3d','sciencedates','histutils','astrometry_azel', 'morecvutils', 'gridaurora', 'lowtran',
]
tests_require=['nose','coveralls']

# %%
import subprocess # need for git
from pathlib import Path
from setuptools import setup,find_packages
      
# leave astropy in here for gaussfitter


#%%
# FIXME: trick is to have them in order from no prereq to full prereq
for p in [
        'https://github.com/scivision/pybashutils',
        'https://github.com/scivision/pyimagevideo',
         'https://github.com/scivision/dmcutils',
          'https://github.com/scivision/themisasi',
        'https://github.com/scivision/transcarread',
        'https://github.com/scivision/gaussfitter',
        'https://github.com/scivision/dascutils',]:

    cwd = Path('..') / p.split('/')[-1]
    print(f'\n {cwd} \n')

    if not cwd.is_dir():
        subprocess.run(['git','clone', p],cwd='..')

    subprocess.call(['git','pull'], cwd=cwd) #in case it was already installed

    subprocess.call(['python','setup.py','develop'], cwd=cwd) #FIXME is there an API to do this?

#%%
setup(name='histfeas',
      packages=find_packages(),
	  description='Feasibility study for HiST auroral tomography system',
	  author='Michael Hirsch, Ph.D.',
	  url='https://github.com/scivision/histfeas',
	  extras_require = {'tiff':['tifffile'],
                        'tests':tests_require},
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',
      ],
      python_requires='>=3.6',
      install_requires=req,
      tests_require=tests_require,
	  )


