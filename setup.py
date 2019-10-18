 #!/usr/bin/env python
import os
import re
import sys
import warnings
import versioneer
#from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

DISTNAME = 'xcape'
LICENSE = 'MIT'
AUTHOR = 'xcape Developers'
AUTHOR_EMAIL = 'rpa@ldeo.columbia.edu'
URL = 'https://github.com/xgcm/xcape'
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: Apache Software License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering',
]

INSTALL_REQUIRES = ['xarray>=0.12.0', 'dask', 'numpy>=1.16']
PYTHON_REQUIRES = '>=3.6'

ext_cape_ml = Extension(name = 'CAPE_CODE_model_lev',
                        sources = ['xcape/CAPE_CODE_model_lev.pyf',
                                   'xcape/CAPE_CODE_model_lev.f90'])
ext_cape_pl = Extension(name = 'CAPE_CODE_pressure_lev',
                        sources = ['xcape/CAPE_CODE_pressure_lev.pyf',
                                   'xcape/CAPE_CODE_pressure_lev.f90'])
ext_bunkers_ml = Extension(name = 'Bunkers_model_lev',
                        sources = ['xcape/Bunkers_model_lev.pyf',
                                   'xcape/Bunkers_model_lev.f90'])
ext_bunkers_pl = Extension(name = 'Bunkers_pressure_lev',
                        sources = ['xcape/Bunkers_pressure_lev.pyf',
                                   'xcape/Bunkers_pressure_lev.f90'])
ext_srh_ml = Extension(name = 'SREH_model_lev',
                        sources = ['xcape/SREH_model_lev.pyf',
                                   'xcape/SREH_model_lev.f90'])
ext_srh_pl = Extension(name = 'SREH_pressure_lev',
                        sources = ['xcape/SREH_pressure_lev.pyf',
                                   'xcape/SREH_pressure_lev.f90'])
ext_stdh_ml = Extension(name = 'stdheight_2D_model_lev',
                        sources = ['xcape/stdheight_2D_model_lev.pyf',
                                   'xcape/stdheight_2D_model_lev.f90'])
ext_stdh_pl = Extension(name = 'stdheight_2D_pressure_lev',
                        sources = ['xcape/stdheight_2D_pressure_lev.pyf',
                                   'xcape/stdheight_2D_pressure_lev.f90'])

DESCRIPTION = "Fast convective parameters for numpy, dask, and xarray"
def readme():
    with open('README.rst') as f:
        return f.read()

setup(name=DISTNAME,
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      license=LICENSE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      classifiers=CLASSIFIERS,
      description=DESCRIPTION,
      long_description=readme(),
      #install_requires=INSTALL_REQUIRES,
      #python_requires=PYTHON_REQUIRES,
      url=URL,
      #packages=['xcape'],
      ext_package='xcape',
      ext_modules = [ext_cape_ml, ext_cape_pl, 
                     ext_bunkers_ml, ext_bunkers_pl, 
                     ext_srh_ml, ext_srh_pl, 
                     ext_stdh_ml, ext_stdh_pl]
      )
