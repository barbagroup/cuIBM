"""
Sets up the version.
"""

import os


_version_major = 0
_version_minor = 3
_version_micro = ''
_version_extra = ''

# construct full version string
_ver = [_version_major, _version_minor]
if _version_micro:
  _ver.append(_version_micro)
if _version_extra:
  _ver.append(_version_extra)
__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'License :: OSI Approved :: MIT License',
               'Operating System :: Unix',
               'Programming Language :: Python']

description = 'snake: post-processing tools for the flying-snake simulations'
# Long description will go up on the pypi page
long_description = """
Snake
=====
Snake is a collection of Python modules used to post-process the numerical
solution of flying-snake simulations using one of the following software:

  * [cuIBM](https://github.com/barbagroup/cuIBM):
    a GPU-based immersed boundary method code;
  * [PetIBM](https://github.com/barbagroup/PetIBM):
    a parallel immersed boundary method code;
  * [IBAMR](https://github.com/IBAMR/IBAMR):
    an adaptive and distributed-memory parallel implementation of the immersed
    boundary (IB) method;
  * IcoFOAM: the incompressible laminar solver of
    [OpenFOAM](http://www.openfoam.org/).

License
=======
``snake`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.
All trademarks referenced herein are property of their respective holders.
Copyright (c) 2016--, Olivier Mesnard, The George Washington University.
"""

NAME = 'snake'
MAINTAINER = 'Olivier Mesnard'
MAINTAINER_EMAIL = 'mesnardo@gwu.edu'
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = 'https://github.com/mesnardo/snake'
DOWNLOAD_URL = ''
LICENSE = 'MIT'
AUTHOR = 'Olivier Mesnard'
AUTHOR_EMAIL = 'mesnardo@gwu.edu'
PLATFORMS = 'Unix'
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['snake',
            'snake.cuibm',
            'snake.petibm',
            'snake.openfoam',
            'snake.ibamr',
            'snake.openfoam',
            'snake.solutions',
            'snake.tests']
PACKAGE_DATA = {'snake': [os.path.join('styles', '*')]}
REQUIRES = ['numpy', 'matplotlib', 'scipy', 'pandas']
