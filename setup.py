#!/usr/bin/env python
from setuptools import setup, find_packages

install_requires = ['numpy', 'h5py', 'scipy', 'pandas', 'xarray', 'matplotlib', 'seaborn', 'astropy',
                    'Wand', 'pathvalidate', 'geopy', 'simplekml',
                    'pymap3d', 'sciencedates', 'histutils', 'astrometry_azel', 'morecvutils', 'gridaurora', 'lowtran',
                    'dascutils', 'transcarread', 'themisasi', 'pyimagevideo',
                    'pybashutils', 'dmcutils',
                    ]
tests_require = ['pytest', 'nose', 'coveralls']


# %%
setup(name='histfeas',
      packages=find_packages(),
      description='Feasibility study for HiST auroral tomography system',
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      version='0.6.0',
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/histfeas',
          extras_require={'tiff': ['tifffile'],
                          'maps': ['cartopy'],
                          'fit': ['gaussfitter'],
                          'tests': tests_require},
      classifiers=[
          'Intended Audience :: Science/Research',
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
          'Programming Language :: Python :: 3.6',
      ],
      python_requires='>=3.6',
      install_requires=install_requires,
      tests_require=tests_require,
      scripts=['EnergyASI.py', 'FigureMaker.py', 'RunHistfeas.py'],
      )
