#!/usr/bin/env python
import os,sys,subprocess
from setuptools import setup

exepath = os.path.dirname(sys.executable)
try:
    subprocess.run([os.path.join(exepath,'conda'),'install','--yes','--file','requirements.txt'])
except Exception as e:
    print('tried conda in {}, but you will need to install packages in requirements.txt  {}'.format(exepath,e))


with open('README.rst','r') as f:
	long_description = f.read()

setup(name='histfeas',
      version='0.1',
	  description='Feasibility study for HiST auroral tomography system',
	  long_description=long_description,
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/histfeas',
	  install_requires=['histutils','lowtran','pymap3d','astrometry_azel','cvutils','gridaurora','transcarread','pyimagevideo','dascutils','pybashutils',
                        'gaussfitter',
                        'tifffile','Wand',
		        'geopy'],
      dependency_links = [
        'https://github.com/scienceopen/histutils/tarball/master#egg=histutils',
        'https://github.com/scienceopen/lowtran/tarball/master#egg=lowtran',
        'https://github.com/scienceopen/pymap3d/tarball/master#egg=pymap3d',
        'https://github.com/scienceopen/astrometry_azel/tarball/master#egg=astrometry_azel',
        'https://github.com/scienceopen/CVutils/tarball/master#egg=CVutils',
        'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora',
        'https://github.com/scienceopen/transcarread/tarball/master#egg=transcarread',
        'https://github.com/scienceopen/gaussfitter/tarball/master#egg=gaussfitter',
        'https://github.com/scienceopen/pyimagevideo/tarball/master#egg=pyimagevideo',
        'https://github.com/scienceopen/dascutils/tarball/master#egg=dascutils',
        'https://github.com/scienceopen/pybashutils/tarball/master#egg=pybashutils'
                          ],
      packages=['histfeas'],
	  )

