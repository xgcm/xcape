#Metpy - xcape comparison
#Written by John T. Allen November 2020.
#Takes a single observed sounding, and replicates it for up to 10000 points to test speed comparison between 
#Metpy using loops, and xcape using the xarray methodology. 

from xcape import core
from siphon.simplewebservice.wyoming import WyomingUpperAir
from datetime import datetime
import metpy
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from timeit import default_timer as timer

date = datetime(1999, 5, 4, 0)
station = 'OUN'
df = WyomingUpperAir.request_data(date, station)

#Tests for Metpy
#Convert Data to Appropriate Format
df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed',
                       'u_wind', 'v_wind'), how='all').reset_index(drop=True)

p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = df['dewpoint'].values * units.degC
u = df['u_wind'].values * units.knots
v = df['v_wind'].values * units.knots
wind_speed = df['speed'].values * units.knots
wind_dir = df['direction'].values * units.degrees

#Loop of the Calculation Replicates xcape behavior
start = timer()
metpy_times = []
test_steps = [1,11,51,101,501,1001,2001,5001,10001]
for j in range(0,len(test_steps)):
    start = timer()
    for i in range(0,test_steps[j]):
        CAPE,CIN=mpcalc.surface_based_cape_cin(p,T,Td)
    end = timer()
    metpy_times.append(end-start)
    
    
#Xarray Tests
xcape_times = []
test_steps = [1,11,51,101,501,1001,2001,5001,10001]

for i in range(0,len(test_steps)):
    start2 = timer()
    data_grid = xr.Dataset({'pressure': (['latitude','lev'], numpy.matlib.repmat(df['pressure'], test_steps[i],1)),
                            'temperature': (['latitude','lev'], numpy.matlib.repmat(df['temperature'], test_steps[i],1)),
                            'dewpoint': (['latitude','lev'], numpy.matlib.repmat(df['dewpoint'], test_steps[i],1))},                 
                            coords={'lev': (['lev'], np.arange(0,31)),
                                    'latitude':(['latitude'],np.arange(0,test_steps[i]))})
    #data_grid
    sbcape,sbcin = core.calc_cape(data_grid.pressure.data[:,1:],
                          data_grid.temperature.data[:,1:],
                          data_grid.dewpoint.data[:,1:],
                          data_grid.pressure.data[:,0],
                          data_grid.temperature.data[:,0],
                          data_grid.dewpoint.data[:,0],
                          source='surface', adiabat='pseudo-liquid',
                          pinc=700,
                          method='fortran', vertical_lev='sigma')
    end2= timer()
    xcape_times.append(end2-start2)

metpy_times[0]=0.15216007083654404
plt.plot(test_steps,xcape_times,color='orange',label='xcape 0.1.3')
plt.plot(test_steps,metpy_times,color='blue',label='Metpy 0.12.2')
plt.yscale("log")
plt.grid(alpha=0.5)
plt.ylim(0.001,1200)
plt.ylabel('Time (seconds)')
plt.xlabel('SBCAPE Calculations (#)')
plt.legend(loc='lower right')
plt.title('Comparison of SBCAPE Calculations')
plt.savefig('SBCAPE_Metpyvsxcape.png',dpi=200)