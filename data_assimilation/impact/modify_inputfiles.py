import xarray as xr
import sys, os 
import numpy as np
import glob

def modify_inputfiles(f, pm, pn, ntimes, nHIS, spherical, Vtransform, Vstretching, theta_s, theta_b, Tcline, mask_rho, mask_u, mask_v):
    path = '... path/to/files/' + f.split('/')[-1]

    if os.path.isfile(path):
        return path
    else:
        ds = xr.open_dataset(f)

        ds['pm'] = pm
        ds['pn'] = pn
        ds['ntimes'] = ntimes 
        ds['nHIS'] = nHIS
        ds['spherical'] = spherical
        ds['Vtransform'] = Vtransform
        ds['Vstretching'] = Vstretching
        ds['theta_s'] = theta_s
        ds['theta_b'] = theta_b
        ds['Tcline'] = Tcline
        
        if isinstance(mask_rho, xr.DataArray):
            ds['mask_rho'] = mask_rho
            ds['mask_u'] = mask_u
            ds['mask_v'] = mask_v
        else:
            ds['mask_rho'] = xr.DataArray(mask_rho, dims = ('eta_rho', 'xi_rho'))
            ds['mask_u'] = xr.DataArray(mask_u, dims = ('eta_u', 'xi_u'))
            ds['mask_v'] = xr.DataArray(mask_v, dims = ('eta_v', 'xi_v'))
        
        ds.to_netcdf(path)
        return path


tmpfile = glob.glob('/path/to/file/with/pm_and_pn_variables/norshelf-ana_qck_*')[0]
tmpds = xr.open_dataset(tmpfile)
pm = tmpds.pm
pn = tmpds.pn
tmpds.close()

files = sorted(glob.glob('/path/to/bck/files/we/calculate/impact/from/norshelf-ana_fwd_*_outer0.nc'))

for f in files:
    print(f)
    d = f.split('fwd_')[-1].split('T')[0]
    tmpfile = '/path/to/daifiles/norshelf-ana_dai_{0}T00Z.nc'.format(d)

    tmpds = xr.open_dataset(tmpfile)
    ntimes = tmpds.ntimes.values
    nHIS = tmpds.nHIS.values
    spherical = tmpds.spherical.values
    Vtransform = tmpds.Vtransform.values
    Vstretching = tmpds.Vstretching.values
    theta_s = tmpds.theta_s.values
    theta_b = tmpds.theta_b.values
    Tcline = tmpds.Tcline.values
    mask_rho = tmpds.mask_rho.values
    mask_u = tmpds.mask_u.values
    mask_v = tmpds.mask_v.values

    tmpds.close()
    modify_inputfiles(f, pm, pn, ntimes, nHIS, spherical, Vtransform, Vstretching, theta_s, theta_b, Tcline, mask_rho, mask_u, mask_v)