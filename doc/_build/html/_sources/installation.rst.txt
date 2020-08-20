
Installation
------------

Requirements
^^^^^^^^^^^^

xcape is compatible with python 3. It requires numpy and, optionally,
xarray.

Installation from Conda Forge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install xcape along with its dependencies is via conda
forge::

    conda install -c conda-forge xcape


Installation from Pip
^^^^^^^^^^^^^^^^^^^^^

An alternative is to use pip::

    pip install xcape

This will install the latest release from
`pypi <https://pypi.python.org/pypi>`_.

Installation from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^

xcape is under active development. To obtain the latest development version,
you may clone the `source repository <https://github.com/xgcm/xcape>`_
and install it::

    git clone https://github.com/xgcm/xcape.git
    cd xcape
    python setup.py install

or simply::

    pip install git+https://github.com/xgcm/xcape.git

Users are encouraged to `fork <https://help.github.com/articles/fork-a-repo/>`_
xcape and submit issues_ and `pull requests`_.

.. _dask: http://dask.pydata.org
.. _xarray: http://xarray.pydata.org
.. _issues: https://github.com/xgcm/xcape/issues
.. _`pull requests`: https://github.com/xgcm/xcape/pulls
