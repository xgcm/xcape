"""
Metadata for xarray outputs.
"""

def cape_metadata(kwargs):
    """"Return attribute dict for cape xarray output."""

    surf_suffix = '_wrt_surface' if kwargs['source'] == 'surface' else ''
    surf_suffix_long = ' w.r.t. Surface' if kwargs['source'] == 'surface' else ''

    attrs = {
        'CAPE': {
            'standard_name': 'atmosphere_convective_available_potential_energy' + surf_suffix,
            'long_name': 'Convective Available Potential Energy' + surf_suffix_long,
            'units': 'J kg-1',
            'cape_parcel_source': kwargs['source'],
            # 'cape_pressure_increment': kwargs['pinc'],
            # 'cape_vertical_lev': kwargs['vertical_lev'],
            # 'cape_adiabat': kwargs['adiabat']
        },
        'CIN': {
            'standard_name': 'atmosphere_convective_inhibition' + surf_suffix,
            'long_name': 'Convective Inhibition' + surf_suffix_long,
            'units': 'J kg-1',
            'cape_parcel_source': kwargs['source'],
            # 'cape_pressure_increment': kwargs['pinc'],
            # 'cape_vertical_lev': kwargs['vertical_lev'],
            # 'cape_adiabat': kwargs['adiabat']
        },
        'mulev': {
            'long_name': 'Most unstable level index'
        },
        'z_mulev': {
            'long_name': 'Height of most unstable level',
            'units': 'm'
        }
    }

    if kwargs['source'] == 'mixed-layer':
        attrs['CAPE']['mixed_layer_depth'] = f"{kwargs['ml_depth']} m"
        attrs['CIN']['mixed_layer_depth'] = f"{kwargs['ml_depth']} m"

    return attrs
