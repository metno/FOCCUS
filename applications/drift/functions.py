"""
Functions for extracting model variables along a predefined OpenDrift trajectory.

value_along_trajectory() (slower) and model_variables_along_trajectory() (faster) should yield the same results.

Author: Mateusz Matuszak
"""
import numpy as np

def _distance(lat1, lon1, lat2, lon2):
    from math import cos, asin
    "Haversine formula"
    p = np.pi/180
    hav = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p)*cos(lat2*p) * (1-cos((lon2-lon1)*p)) / 2
    return 12742 * asin(np.sqrt(hav))

def _closest(data, target_lat, target_lon):
    return min(data, key=lambda p: _distance(target_lat, target_lon, p[0], p[1]))

def value_along_trajectory(trajectory, DataArray, variable='temperature'):
    from joblib import Parallel, delayed
    from joblib_progress import joblib_progress
    values = []
    
    def single_trajectory(traj):
        i = 0
        val = np.zeros(len(trajectory.time.values))
        for time, lon, lat, z in zip(trajectory.time.values, trajectory.lon.values[traj], trajectory.lat.values[traj], trajectory.z.values[traj]):
            tmp = DataArray.sel(time=time, depth=z, method='nearest')
            tmp = tmp.where((DataArray.lon > lon - 0.05) & 
                            (DataArray.lon < lon + 0.05) & 
                            (DataArray.lat > lat - 0.025) & 
                            (DataArray.lat < lat + 0.025), drop=True)
            nearest_grid = _closest([[ll, lt] for ll, lt in zip(tmp.lat.values.flatten(), tmp.lon.values.flatten())], lat, lon)
            nearest_lon, nearest_lat = nearest_grid[1], nearest_grid[0]
            tmp = tmp.where((tmp.lon == nearest_lon) & (tmp.lat == nearest_lat), drop=True)
            val[i] = tmp[variable].values[0,0]
            i+=1
        return val
    with joblib_progress("Finding values along trajectory...", total=(len(trajectory.trajectory))):
        values.append(Parallel(n_jobs=8)(delayed(single_trajectory)(traj) for traj in range(len(trajectory.trajectory))))
    values = np.array(values[0])
    return values

def model_variables_along_trajectory(traj='sample_file.nc', model='https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be', variable='temperature', trajectory_nr=0, model_res=800, constant_time=None, times=None):
    import pyproj
    import pyresample as py
    import xarray as xr
    traj = xr.open_dataset(traj)
    ds = xr.open_dataset(model)
    crs = ds.projection_stere.proj4
    if variable == 'abs_vel':
        u, v = ds['u_eastward'], ds['v_northward']
        ds = np.sqrt(u**2+v**2)
    else:
        ds = ds[variable]

    ds_geo = py.geometry.GridDefinition(lons=ds.lon, lats=ds.lat)
    traj_geo = py.geometry.SwathDefinition(lons=traj.lon, lats=traj.lat, crs=crs)

    _, valid_output_index, index_array, distance_array = \
                                py.kd_tree.get_neighbour_info(
                                    source_geo_def=ds_geo,
                                    target_geo_def=traj_geo,
                                    radius_of_influence=model_res,
                                    neighbours=1)
    
    index_array_2d = np.unravel_index(index_array, ds_geo.shape)
    if times is None:
        iteration_times = traj.time.values
    else:
        iteration_times = times
    values_array = np.zeros(len(iteration_times))
    i = 0
    for time, z in zip(iteration_times, traj.z.values[trajectory_nr]):
        X, Y = index_array_2d[1][i], index_array_2d[0][i]
        if constant_time is not None:
            value = ds.sel(time=constant_time)
        else:
            value = ds.sel(time=time)
        value = value.sel(depth=z, method='nearest')
        value = value.isel(X=X, Y=Y)
        values_array[i] = value.values
        i+=1

    return values_array



if __name__ == '__main__':
    model_variables_along_trajectory()