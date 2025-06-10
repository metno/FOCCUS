from subprocess import call
from netCDF4 import Dataset
import numpy as np


def create_adsenfile(fwdfile, adsendir):
    # Split file name in directory path and name
    adsfile = '/'.join([adsendir, fwdfile.split('/')[-1].replace('fwd','ads')])#.split('-')[0]
    print(adsfile) 
    cmdstring = 'ncgen -k4 -o ' + adsfile + ' -b adsen.cdl'

    call(cmdstring, shell=True)
    
    fout = Dataset(adsfile,'r+',)
    f = Dataset(fwdfile,'r',)
    N = len(f.dimensions['ocean_time'])
    # Copy time info
    print('Copying time')
    fout.variables['ocean_time'][:] = f.variables['ocean_time'][:]
    

    # Copy stretching info
    copy_stretch = ['spherical', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'Tcline', 'hc', 's_rho', 's_w', 'Cs_r', 'Cs_w']
    for var in copy_stretch:
        print('Copying variable: ' + var)
        fout.variables[var][:] = f.variables[var][:]  
        
    # Copy grid info
    copy_grid = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']
    for var in copy_grid:
        print('Copying variable: ' + var)
        fout.variables[var][:] = f.variables[var][:]
        
        
    # Initialize scope from masks
    print('Setting scope for rho/u/v')
    fout.variables['scope_rho'][:,:] = f.variables['mask_rho'][:,:]
    fout.variables['scope_u'][:,:] = f.variables['mask_u'][:,:]
    fout.variables['scope_v'][:,:] = f.variables['mask_v'][:,:]
    
    # Set 2D state variables to zero
    copy_var2D = ['zeta', 'ubar', 'vbar']
    for var in copy_var2D:
        print('Initializing variable: ' + var)
        fout.variables[var][:] = 0
        
    # Set 3D state variables to zero
    copy_var3D = ['u', 'v', 'temp', 'salt']
    for var in copy_var3D:
        print('Initializing variable: ' + var)
        fout.variables[var][:] = 0
    
    f.close()
    fout.close()
    
    return adsfile
