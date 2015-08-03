#!/usr/bin/env python3

from setuptools import setup 

with open('README.rst') as f:
	long_description = f.read()
	
setup(name='histfeas',
      version='0.1',
	  description='Feasibility study for HiST auroral tomography system',
	  long_description=long_description,
	  author='Michael Hirsch',
	  author_email='hirsch617@gmail.com',
	  url='https://github.com/scienceopen/histfeas',
	  install_requires=['histutils','lowtran','pymap3d','astrometry_azel','CVutils','gridaurora','transcarread',
                        'tifffile','seaborn',
                        'numpy','pytz','pandas','scipy','astropy'],
      dependency_links = ['https://github.com/scienceopen/histutils/tarball/master#egg=histutils',
                          'https://github.com/scienceopen/lowtran/tarball/master#egg=lowtran',
                          'https://github.com/scienceopen/pymap3d/tarball/master#egg=pymap3d',
                          'https://github.com/scienceopen/astrometry_azel/tarball/master#egg=astrometry_azel',
                          'https://github.com/scienceopen/CVutils/tarball/master#egg=CVutils',
                          'https://github.com/scienceopen/gridaurora/tarball/master#egg=gridaurora',
                          'https://github.com/scienceopen/transcarread/tarball/master#egg=transcarread'
                          ],
      packages=['histfeas'],
	  )
