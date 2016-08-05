#!/usr/bin/env python
import subprocess
from setuptools import setup

try:
    subprocess.call(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    pass

setup(name='histfeas',
      packages=['histfeas'],
	  description='Feasibility study for HiST auroral tomography system',
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/histfeas',
	  install_requires=['themisasi','histutils','lowtran','pymap3d','astrometry_azel',
			'cvutils','gridaurora','transcarread','pyimagevideo','dascutils','pythonutils',
                        'gaussfitter',
                        'Wand','pathvalidate',
		         'geopy','simplekml'],
	  extras_require = {'tifffile':['tifffile']},
      dependency_links = [
        'https://github.com/scienceopen/themisasi/tarball/master#egg=themisasi',
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
        'https://github.com/scienceopen/pybashutils/tarball/master#egg=pythonutils'
                          ],
	  )

