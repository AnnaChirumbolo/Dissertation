import numpy as np
import sys
import xarray as xr #xarray to read all types of formats
from affine import Affine
import glob 

# create 1d array of lat/lons

llats = np.arange(89.75,-90.25,-.5); llons = np.arange(-179.75,180.25,.5)

ny = llats.size
nx = llons.size

coords = {'lat': (['lat'],llats,{'units':'degrees_north','long_name':'latitude'}),
          'lon': (['lon'],llons,{'units':'degrees_east','long_name':'longitude'})}

attrs = []
attrs.append({'_FillValue':-9999.,'units':'mm2.kg-1','long_name':'specific leaf area'})
attrs.append({'_FillValue':-9999.,'units':'mg.g-1','long_name':'leaf nitrogen'})
attrs.append({'_FillValue':-9999.,'units':'mg.g-1','long_name':'leaf phosphorus'})

var_names = ['sla','lnm','lpm']
std_names = ['sla_std','lnm_std','lpm_std']



# create 2d grid of latitudes/longitudes
lon2d,lat2d = np.meshgrid(llons,llats)

# create template xarray
data_vars={}
for ff in range(0,3):
	data_vars[var_names[ff]] = (['lat','lon'],np.zeros([ny,nx])*np.nan,attrs[ff])
	data_vars[std_names[ff]] = (['lat','lon'],np.zeros([ny,nx])*np.nan,attrs[ff])
ds_new = xr.Dataset(data_vars=data_vars,coords=coords)


sequence = ('sla.csv','lnm.csv','lpm.csv')
for ff, filename in enumerate(sequence, start=0):
	print('filename')
	ds = np.loadtxt(filename,delimiter=",", skiprows=1)
	var_2d = np.empty((len(llats),len(llons),))
	var_2d[:] = np.nan

	var_std_2d = np.empty((len(llats),len(llons),))
	var_std_2d[:] = np.nan

	lat = ds[:,1]
	lon = ds[:,0]
	var = ds[:,2]
	var_std = ds[:,3]

	i = 0
	while i < len(ds):
		#print("lon = "+str(ds[:,0][i])+", lat = "+str(ds[:,1][i])+", sla = "+str(ds[:,2][i])+", sla_sd = "+str(ds[:,3][i]))
		if i%10000 == 0:
			print ('\r\tprocessing row %i of %i' % (i,len(ds)))
		# regrid
		slc = (np.abs(lat[i]-lat2d)<.25)*(np.abs(lon[i]-lon2d)<.25)
		lat_tmp, lon_tmp = np.where(slc==1)
	
		#print lat_tmp, lon_tmp
		#print sla_2d.shape
		#print sla.shape
	
		var_2d[lat_tmp[0], lon_tmp[0]]= var[i]
		var_std_2d[lat_tmp[0], lon_tmp[0]]= var_std[i]
		
		i += 1

	ds_new[var_names[ff]].values = var_2d.copy()
	ds_new[std_names[ff]].values = var_std_2d.copy()


nc_file = 'Butler_Leaftraits_Processed.nc'
ds_new.to_netcdf(path=nc_file)

#used cdo to convert to desired spacial resolution (0.25, 0.5 and 1 degree)
