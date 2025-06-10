# This script generates a monthly z-depth file

import xarray as xr
import numpy as np
import xroms
from utils import get_boundary_slice
#ignore some xroms warnings
import warnings
warnings.filterwarnings("ignore")

#-------------------------------------------------------------------------
dir_avr='./' #SPECIFY OUTPUT DIR HERE
outer = 0

drop_vars = ['AKs','AKt','AKv','DU_avg1','DU_avg2','DV_avg1','DV_avg2',
             'u_obc','ubar','ubar_obc','v_obc','vbar','vbar_obc','zeta_obc',
             'lat_u','lat_v','lon_u','lon_v']
# TODO: should not drop the remaining _obc variables

#-------------------------------------------------------------------------
def get_ds(dir,file):
    ds = xr.open_dataset(dir + file,drop_variables=['ocean_time']+drop_vars)
    # NOTE: ocean_time is shifted 1d back in time, so use 'time' and rename it to 'ocean_time'
    # (this is because of the post-processing using nco)
    # Remove the ocean_time dimension if it exists and is not needed
    if 'ocean_time' in ds.dims and 'ocean_time' not in ds.variables:
        ds = ds.reset_coords('ocean_time', drop=True)
    # Now rename 'time' to 'ocean_time'
    ds = ds.rename({'time': 'ocean_time'})
    return ds

#-------------------------------------------------------------------------
# Get the monthly dataset & extract the obc variables & reshape them
ds_list = []
for i in range(1,13):
    # Get the file name
    file = f'norshelf-ana_fwd_avr_2024{i:02d}_outer{outer}.nc'
    ds_list.append(get_ds(dir_avr,file))
ds = xr.concat(ds_list, dim='ocean_time')
print("Got the monthly dataset")

# get obc variables
obc_vars = []
for var in ds.data_vars:
    if 'obc' in var:
        obc_vars.append(var)
ds_obc = ds[obc_vars].copy()

# Then remove the obc variables from the original dataset
ds = ds.drop_vars(obc_vars)

#------------
# Reshape obc variables into 4D arrays which then can be merged with ds
# _obc variables will not be converted by xroms since their dims/axes are strange

# initialise a ds with same dimensions as ds
ds_obc_XY = ds.copy().drop_vars(ds.data_vars)
# TODO: to save memory, this could be two datasets each for X and Y dims
# which have only two indexes (X: east-west, Y: south-north)
# This would avoid saving all the extra nans...

# Get an array of the ocean_time dimension
ntime = ds_obc.ocean_time.values.shape[0]
itime = np.arange(ntime)

for var_name in ds_obc.data_vars:
    # add the all dimensions to the variable (fill with nans in the center of domain)
    #('ocean_time', 's_rho', 'eta_rho', 'xi_rho')=(12, 42, 351, 901)
    var_full = np.zeros((ntime,42,351,901))
    var_name = var_name.split('_')[0]

    for dir,idir in zip(['east','west'],[0,-1]):
        # NOTE: var variables such as u, v, the length of the indexes are different...see utils.py. Fix later
        var_slice = get_boundary_slice(ds_obc,var_name,direction=dir,itime=itime,obc=True,obc_adjust=0)
        var_full[:,:,idir,:] = var_slice
        
    for dir,idir in zip(['south','north'],[0,-1]):
        var_slice = get_boundary_slice(ds_obc,var_name,direction=dir,itime=itime,obc=True,obc_adjust=0)
        var_full[:,:,:,idir] = var_slice[:,:,:351]
        
    ds_obc_XY[var_name+'_obc'] = (('ocean_time','s_rho','eta_rho', 'xi_rho'),var_full)

# set ocean_time to the same as ds
ds_obc_XY['ocean_time'] = ds.ocean_time
print("Reshaped the obc variables")

# Done, ds_obc_XY now contains only values on the boundary, nan elsewhere
del ds_obc

# Then merge the obc variables into ds
ds = xr.merge([ds, ds_obc_XY])
del ds_obc_XY

#-------------------------------------------------------------------------
# Use xroms to get the z-coordinates

# Need this file to get the 'pm' and 'pn' variables since they have been removed 
# in the processing of the fwd-files
tmpfile = dir_avr + 'norshelf-ana_qck_20250518T00Z_000.nc'
tmpds = xroms.open_netcdf(tmpfile, chunks={"ocean_time": 1}, Vtransform=2)
pm = tmpds.pm
pn = tmpds.pn
tmpds.close()
# delete tmpds
del tmpds

# Get the z-coordinates
# Dropping variables that xroms produces that are not needed
var_drop_list = ['Cs_r', 'Cs_w', 'hc', 'pm', 'pn', 'spherical', '3d', 'dx', 'dx_u', 'dx_v', 'dx_psi', 'dy', 'dy_u', 'dy_v', 'dy_psi', 'dz', 'dz_w', 'dz_u', 'dz_w_u', 'dz_v', 'dz_w_v', 'dz_psi', 'dz_w_psi', 'dA', 'rho0']
print("Starting to get the z-coordinates with xroms")
def get_dz(ds, pm, pn, ntime=None):
    if ntime is None:
        print('No ntime given, using length of ocean_time')
        ntime = ds.dims['ocean_time']
    dz_list = []
    for i in range(ntime):
        print(f'Getting z-coordinates for time step {i+1}/{ntime}')
        # Select the time step
        if 'ocean_time' in ds.dims:
            ds_it = ds.isel(ocean_time=i)
        else:
            ds_it = ds
        ds_it.salt.shape
        # add pm and pn
        ds_it['pm'] = pm
        ds_it['pn'] = pn
        # set ocean_time to the first date in each month using the counter i
        new_date = f'2024-{i+1:02d}-01T00:00:00' 
        ds_it['ocean_time'] = np.datetime64(new_date)
        print(ds_it.ocean_time.values)

        # Use xroms to get the z-grid
        dzi, grid = xroms.roms_dataset(ds_it, Vtransform=2)

        # append
        dz_list.append(dzi.drop_vars(var_drop_list, errors='ignore'))
    print("Concatenating...")
    dz = xr.concat(dz_list, dim='ocean_time')

    return dz, grid

dz, grid0 = get_dz(ds, pm, pn)
print(dz.ocean_time.values)

# Set ocean_time to the same as ds

#--------------------------------------------------------------------------
# save dz to a file
outfile = dir_avr + f'norshelf-ana_fwd_avr_zdepth_202401-12_outer{outer}.nc'
print("Saving the z-coordinates to file "+outfile)
dz.to_netcdf(outfile, mode='w', format='netcdf4', encoding={'ocean_time': {'units': 'seconds since 1970-01-01 00:00:00', 'calendar': 'proleptic_gregorian'}})
print("DONE")
