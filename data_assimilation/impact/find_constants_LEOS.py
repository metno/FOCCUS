import numpy as np
import sys
import seawater as sw
import os, glob
import xroms
from datetime import datetime
import xarray as xr
from typing import Sequence
from latlon import ll2xy #Python script
import warnings
warnings.filterwarnings("ignore")

def section(A: xr.Dataset, X: Sequence[float], Y: Sequence[float]) -> xr.Dataset:
    """Slice a Dataset along a vertical section"""

    # Make temporary DataArrays for initial interpolations
    X0 = xr.DataArray(X, dims=["distance"])
    Y0 = xr.DataArray(Y, dims=["distance"])
    # Distance between section points
    pm = A.pm.interp(xi_rho=X0, eta_rho=Y0).values
    pn = A.pn.interp(xi_rho=X0, eta_rho=Y0).values
    dX = 2 * np.diff(X) / (pm[:-1] + pm[1:])
    dY = 2 * np.diff(Y) / (pn[:-1] + pn[1:])
    dS = np.sqrt(dX * dX + dY * dY)
    # Cumulative distance along the section
    distance = np.concatenate(([0], np.cumsum(dS))) / 1000.0  # unit = km
    X0["distance"] = distance
    Y0["distance"] = distance
    # Interpolate to the section making an intermediate Dataset
    B0 = A.interp(xi_rho=X0, eta_rho=Y0)

    # Initialize the proper section Dataset
    B = xr.Dataset(dict(xi_rho=X, eta_rho=Y, s_rho=A.s_rho))

    # Weights for trapezoidal integration
    V = 0.5 * (np.concatenate(([0], dS)) + np.concatenate((dS, [0])))
    B["dS"] = xr.DataArray(V, dims=["distance"], coords=dict(distance=B0["distance"].values))

    # Scalar data variables
    for var in ["h", "temp", "salt", "zeta", 'sig0', 'mld']:
        print(var)
        if var in B0:
            B[var] = B0[var]

    # Velocity
    if "u" in A:
        B["u"] = A.u.interp(xi_u=X0, eta_rho=Y0)
        B = B.drop(["xi_u"])
    if "v" in A:
        B["v"] = A.v.interp(xi_rho=X0, eta_v=Y0)
        B = B.drop(["eta_v"])
        
    # Normal velocity ++
    if "u" in A:
        dX = diff2(X)
        dY = diff2(Y)
        norm = np.sqrt(dX * dX + dY * dY)
        # Unit normal vector
        nX = -dY / norm
        nY = dX / norm
        B["u_norm"] = B["u"] * nX + B["v"] * nY

    # More vertical structure
    print(B.coords)
    #B.coords["z_w0"] = B0.z_w0
    V = B0.z_w0.values[1:, :] - B0.z_w0.values[:-1, :]
    #B["dZ"] = xr.DataArray(V, dims=["s_rho", "distance"])
    print(V.shape)
    print(B["temp"].values.shape)
    print(B0["temp"].values.shape)
    #B["dZ"] = xr.DataArray(V, dims=["s_rho", "distance"], coords=dict(distance=B["distance"].values, s_rho=B["s_rho"].values))
    B["dZ"] = B["temp"]*5 
    B["dZ"].values[:] = V
    B["area"] = B.dZ * B.dS  # Unit = m**2

    return B

def diff2(X: np.ndarray) -> np.ndarray:
    """Central differences of a sequence"""
    n = len(X)
    Y = np.empty(n + 2, dtype=X.dtype)
    Y[1:-1] = X
    Y[0] = X[0]
    Y[-1] = X[-1]
    return Y[2:] - Y[:-2]


def linearEOS_constants(ds, X, Y, saveto, file, D=False):
    """ Input dataset must be opened by xroms and cannot be as large as the whole grid.
        Set D as the lowest depth you want to calculate constants for. """

    cyclenum = file.split('_fwd_')[-1].split('T')[0]
    dssec = section(ds, X, Y)
    latitude = dssec.lat_rho

    tot_alpha = np.zeros_like(dssec.temp)
    tot_beta = np.zeros_like(tot_alpha)
    tot_T0 = np.zeros_like(tot_alpha)
    tot_S0 = np.zeros_like(tot_alpha)
    tot_dens0 = np.zeros_like(tot_alpha)
        
    N = len(np.shape(tot_alpha))
    print('Shape of input temp/salt:', np.shape(tot_alpha))
    
    temperature = dssec.temp
    salinity = dssec.salt
    density = xroms.potential_density(temperature, salinity, z=0)
    depth = dssec.z_rho0
    
    #for z in range(0, len(ds.s_rho)):
    pressure = sw.eos80.pres(-depth.values, latitude)        
    alpha = sw.eos80.alpha(salinity.values, temperature.values, pressure, pt=True) 
    beta = sw.eos80.beta(salinity.values, temperature.values, pressure, pt=True) 

    tot_alpha[:,:] = alpha
    tot_beta[:,:] = beta
    tot_T0[:,:] = temperature
    tot_S0[:,:] = salinity
    tot_dens0[:,:] = density

    if D:
        tot_alpha[:,:] = np.where(depth >= -D, tot_alpha[:,:], np.nan)
        tot_beta[:,:] = np.where(depth >= -D, tot_beta[:,:], np.nan)
        tot_T0[:,:] = np.where(depth >= -D, tot_T0[:,:], np.nan)
        tot_S0[:,:] = np.where(depth >= -D, tot_S0[:,:], np.nan)
        tot_dens0[:,:] = np.where(depth >= -D, tot_dens0[:,:], np.nan)
    
    if os.path.isfile(saveto):
        with open(saveto, 'a') as f:
            f.write('New cycle: ' + cyclenum + '\n')
            f.write('Alpha:' + str(np.nanmean(tot_alpha)) + '\n')
            f.write('Beta:' + str(np.nanmean(tot_beta)) + '\n')
            f.write('T0:' + str(np.nanmean(tot_T0)) + '\n')
            f.write('S0:' + str(np.nanmean(tot_S0)) + '\n')
            f.write('dens0:' + str(np.nanmean(tot_dens0)) + '\n\n')
    else:
        with open(saveto, 'w') as f:
            f.write('New cycle: ' + cyclenum + '\n')
            f.write('Alpha:' + str(np.nanmean(tot_alpha)) + '\n')
            f.write('Beta:' + str(np.nanmean(tot_beta)) + '\n')
            f.write('T0:' + str(np.nanmean(tot_T0)) + '\n')
            f.write('S0:' + str(np.nanmean(tot_S0)) + '\n')
            f.write('dens0:' + str(np.nanmean(tot_dens0)) + '\n\n')

    return np.array([np.nanmean(tot_alpha), np.nanmean(tot_beta), np.nanmean(tot_T0), np.nanmean(tot_S0), np.nanmean(tot_dens0)])

def linearEOS(constants, T, S): 
    """ Input array as [alpha, beta, T0, S0, rho0] """
    R0 = constants[4] #1027
    T0 = constants[2] #10
    S0 = constants[3] #35
    a = constants[0] #1.7e-4
    b = constants[1] #7.6e-4
    
    return R0 - R0*a*(T-T0) + R0*b*(S-S0)


reftime = datetime(1970, 1, 1)
files = sorted(glob.glob('/lustre/storeB/users/siljeci/FOCCUS/norshelf*'))
print(files)

## saveto = TXT FILE WHERE YOU WANT TO SAVE THE CONSTANTS

latSvinoy = np.array([62.367, 62.485, 62.602, 62.720, 62.780, 62.838, 62.897, 62.957, 
                      63.073, 63.192, 63.310])#, 63.427, 63.663, 63.898, 64.135, 64.370])

lonSvinoy = np.array([5.2, 4.945, 4.690, 4.433, 4.305, 4.175, 4.047, 3.917, 3.657, 3.395,
                      3.133])#, 2.870, 2.338, 1.802, 1.262, 0.728])



cycles = []
if os.path.isfile(saveto):
    with open(saveto, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('New cycle:'):
                cycles.append(int(l.split(': ')[-1]))
print(cycles)

res_list = []
for f in files: 
    num = int(f.split('fwd_')[-1].split('T')[0])
    print(num)

    if num in cycles:
        print('Exists')
        continue
    print('Calculate')

    with xroms.open_netcdf(f, chunks={"ocean_time": 1}, Vtransform=2) as Inp:
        Inp, grid = xroms.roms_dataset(Inp, add_verts=False, include_Z0=True, Vtransform=2)
        
        Inp0 = Inp.mean(dim='ocean_time')
        
        x0, y0 = ll2xy(Inp0, lonSvinoy[0], latSvinoy[0])
        x1, y1 = ll2xy(Inp0, lonSvinoy[-1], latSvinoy[-1])

        i0, j0, i1, j1 = [int(round(v)) for v in [x0, y0, x1, y1]]

        Npoints = 30
        X = np.linspace(x0, x1, Npoints)
        Y = np.linspace(y0, y1, Npoints)
        
        res = linearEOS_constants(Inp0, X, Y, saveto, f, D=98)
        print(res)
        res_list.append(res)
    
