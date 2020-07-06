.. _tutorial:

=============
Tutorial
=============

After installation, **xcape** can be imported in standard python syntax, with the primary functions included in the core element.

  .. code:: python
     >>> import xcape
     or
     >>>from xcape import core

Data can either be model-level or pressure-level based, but will require some structuring using xarray and dask to ensure
that full parallelization of processes is possible. Input data can be parsed from any dataset, using an xr.open_mfdataset
instance and ensuring that parallel and chunk flags are set appropriately for the data being used. An example script for
reading the ERA-Interim reanalysis dataset is provided to illustrate this. Excess variables can be dropped before lazy loading
to improve script efficiency and ensure that RAM bottlenecking isnt a major issue. Required variables are:
1. Near surface temperature, pressure, specific humidity, zonal and meridional wind, and geopotential height if available. 
2. 3D fields of the same parameters. 

Note that if dewpoint temperature is available, this parameter should be also included in the initial data arrays, but otherwise can be
calculated using the included script. It is suggested that the above be of a scripted reading form to ensure that data can be appropriately chunked to suit the system being used for calculation, for example for the provided script:

  .. code:: python
   >>> data, surface = read_reanal(year, month)

This can then be incorporated into a function to calculate convective available potential energy (CAPE):

  .. code:: python 
   >>> def calc_and_save_cape(year, month):
          yearmo = f'{year}{month:02d}'
          data, surface = read_reanal(year, month)
          from xcape import core

          cas,cis = core.calc_cape(data.p.data, data.t.data, data.td.data, 
                  surface.p.data, surface.t.data, surface.td.data,
                  source ='surface',
                  pinc = 700,
                method='fortran',
                vertical_lev='sigma')
          surface['cape_s']=(('time','latitude','longitude'),cas)
          surface['cin_s']=(('time','latitude','longitude'),cis)
   >>> def main():
          calc_and_save_cape(year, month)
   >>> if __name__== "__main__":
          main() 

Note that there are many customizable inputs in the operation calc_cape, including modifying the integration interval for precision,
specifying the vertical level form of the input, and the source of the initial parcel, and whether to use the fortran or numpy implementation (note fortran is preferred for speed). See documentation for further details. 

Note that for large datasets, then an output format of **zarr** is preferable, as it allows for a better parallelization for subsequent
analysis.   
