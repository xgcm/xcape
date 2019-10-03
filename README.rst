xcape: Fast convective parameters for numpy, dask, and xarray
=============================================================

|pypi| |conda forge| |Build Status| |codecov| |docs| |license|

For more information, including installation instructions, read the full
`xcape documentation`_.

.. |conda forge| image:: https://anaconda.org/conda-forge/xcape/badges/version.svg
   :target: https://anaconda.org/conda-forge/xcape
.. |DOI| image:: https://zenodo.org/badge/41581350.svg
   :target: https://zenodo.org/badge/latestdoi/41581350
.. |Build Status| image:: https://travis-ci.org/xgcm/xcape.svg?branch=master
   :target: https://travis-ci.org/xgcm/xcape
   :alt: travis-ci build status
.. |codecov| image:: https://codecov.io/github/xgcm/xcape/coverage.svg?branch=master
   :target: https://codecov.io/github/xgcm/xcape?branch=master
   :alt: code coverage
.. |pypi| image:: https://badge.fury.io/py/xcape.svg
   :target: https://badge.fury.io/py/xcape
   :alt: pypi package
.. |docs| image:: http://readthedocs.org/projects/xcape/badge/?version=latest
   :target: http://xcape.readthedocs.org/en/stable/?badge=latest
   :alt: documentation status
.. |license| image:: https://img.shields.io/github/license/mashape/apistatus.svg
   :target: https://github.com/xgcm/xcape
   :alt: license


Compiling Fortran Code
======================

 f2py  -c ALLCAPELOOP.pyf ALLcapecalcLOOP.f90

