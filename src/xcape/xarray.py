"""
Xarray API for xcap.
"""

import xarray as xr
import numpy as np

from .core import calc_cape as _calc_cape
from .core import calc_srh as _calc_srh


def calc_cape


srh_rm, srh_lm = core.calc_srh(data.p.data, data.t.data, data.td.data, data.u.data, data.v.data,
             surface.p.data, surface.t.data, surface.td.data, surface.u.data, surface.v.data,
                               depth=3000,
                               output_var='srh',
           vertical_lev='sigma')
surface['srh_rm3km']=(('time','lat','lon'),srh_rm)
surface['srh_lm3km']=(('time','lat','lon'),srh_lm)

srh_rm1, srh_lm1 = core.calc_srh(data.p.data, data.t.data, data.td.data, data.u.data, data.v.data,
             surface.p.data, surface.t.data, surface.td.data, surface.u.data, surface.v.data,
                               depth=1000,
                               output_var='srh',
           vertical_lev='sigma')
surface['srh_rm1km']=(('time','lat','lon'),srh_rm1)
surface['srh_lm1km']=(('time','lat','lon'),srh_lm1)


cas,cis = core.calc_cape(data.p.data, data.t.data, data.td.data,
          surface.p.data, surface.t.data, surface.td.data,
          source ='surface',
          pinc = 700,
        method='fortran',
        vertical_lev='sigma')
surface['cape_s']=(('time','lat','lon'),cas)
surface['cin_s']=(('time','lat','lon'),cis)


cas1,cis1 = core.calc_cape(data.p.data, data.t.data, data.td.data,
              surface.p.data, surface.t.data, surface.td.data,
              source ='mixed-layer',
              ml_depth = 500,
              pinc = 700,
            method='fortran',
            vertical_lev='sigma')
surface['cape_ml500']=(('time','latitude','longitude'),cas1)
surface['cin_ml500']=(('time','latitude','longitude'),cis1)

# plt.figure()
# surface.srh_rm.mean('time').plot()
ds_to_save = surface[['cape_s','cin_s',
                      'cape_ml500','cin_ml500',
                      'srh_rm3km', 'srh_lm3km',
                      'srh_rm1km', 'srh_lm1km']]

ds_to_save = ds_to_save.chunk({'time': 150, 'lon': -1, 'lat': -1})
ds_to_save.attrs = data.attrs
ds_to_save.attrs['Generating_code'] = 'https://github.com/xgcm/xcape'
ds_to_save.attrs['Generating_code_version'] = 'v0.1.1'
ds_to_save.attrs['Generating_day'] = date.today().strftime('%d %b %Y')

ds_to_save['cape_s'].attrs = {'long_name': 'Surface Convective Available Potential Energy',
                         'standard_name': 'SBCAPE',
                         'units':'J kg-1'}
ds_to_save['cape_ml500'].attrs = {'long_name': 'Convective Available Potential Energy 500m ml-depth',
                             'standard_name': 'MLCAPE500',
                             'units':'J kg-1'}
ds_to_save['cin_ml500'].attrs = {'long_name': 'Convective Inhibition 500m ml-depth',
                             'standard_name': 'MLCIN500',
                             'units':'J kg-1'}
ds_to_save['cin_s'].attrs = {'long_name': 'Surface Convective Inhibition',
                             'standard_name': 'SBCIN',
                             'units':'J kg-1'}
ds_to_save['srh_rm3km'].attrs = {'long_name': 'Right moving Storm Relative Helicity 0-3km',
                             'standard_name': 'SRH3RM',
                             'units':'m2 s-2'}
ds_to_save['srh_lm3km'].attrs = {'long_name': 'Left moving Storm Relative Helicity 0-3km',
                             'standard_name': 'SRH3LM',
                             'units':'m2 s-2'}
ds_to_save['srh_rm1km'].attrs = {'long_name': 'Right moving Storm Relative Helicity 0-1km',
                             'standard_name': 'SRH1RM',
                             'units':'m2 s-2'}
ds_to_save['srh_lm1km'].attrs = {'long_name': 'Left moving Storm Relative Helicity 0-1km',
                             'standard_name': 'SRH3LM',
                             'units':'m2 s-2'}
