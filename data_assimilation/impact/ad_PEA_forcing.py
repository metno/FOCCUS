"""
This script calculates the adjoint forcing files (using PEA metric) for a given cycle number and outputs the results to a specified directory.
Update directories / file names as needed.
"""

import xarray as xr
import sys, os 
from netCDF4 import Dataset 
import roppy 
import numpy as np
from interpolator import interpolator
from create_adsenfile import create_adsenfile
import subprocess

def mean_rho(constants_file):
    with open(constants_file, 'r') as f:
        for line in f.readlines():
            if not ':' in line:
                continue
            var,val = line.split(':')
            if var == 'Alpha':
                a = float(val)
            elif var == 'Beta':
                b = float(val)
            elif var == 'T0':
                T0 = float(val)
            elif var == 'S0':
                S0 = float(val)
            elif var == 'dens0':
                R0 = float(val)
    return a, b, T0, S0, R0
                
            

def PEA_forcing(cyclenum, fileDir, constants_file, maxdepth ):
    
    # Coordinates of the points we are interested in: 
    lats = np.array( [62.367, 62.39951724, 62.43203448, 62.46455172, 62.49706897, 62.52958621, 62.56210345, 62.59462069, 62.62713793, 62.65965517 ,62.69217241 ,62.72468966 ,62.7572069, 62.78972414, 62.82224138, 62.85475862, 62.88727586, 62.9197931, 62.95231034, 62.98482759, 63.01734483, 63.04986207, 63.08237931, 63.11489655, 63.14741379, 63.17993103, 63.21244828, 63.24496552, 63.27748276, 63.31])

    lons = np.array([5.2, 5.13101715, 5.06203429, 4.99305144, 4.92344291, 4.85403065, 4.78450732, 4.71461734, 4.64453494, 4.57431438, 4.50364706, 4.43342937, 4.36277103, 4.2916826, 4.22068756, 4.14938761, 4.07809246, 4.00632379, 3.93487649, 3.8628644, 3.79040918, 3.71786269, 3.6453162, 3.57276972, 3.50018664, 3.42703098, 3.35387186, 3.28024791, 3.20663875, 3.13301153])

    # The pre-calculated constants:
    a, b, T0, S0, R0 = mean_rho(constants_file)
    print('The constants (R0, a, b): {}, {}, {}'.format(R0,a,b))
    Inp = '/'.join([fileDir, '/norshelf-ana_fwd_{0}T00Z_outer0.nc'.format(cyclenum)])
    Ginp = roppy.SGrid(Dataset(Inp))

    with xr.open_dataset(Inp) as ds:
        dtfac = ds.nHIS.values/ds.ntimes.values

        # Construct output matrices:
        ad_temp = np.zeros_like(ds.temp.values)
        ad_salt = np.zeros_like(ds.salt.values)

        temp_section1 = interpolator(Ginp, ds.isel(ocean_time = 0), 'temp', lons, lats, ad_F=None)
        dis = temp_section1['dis'][0,-1] * 1000
        ad_ind = 9.81 / (maxdepth*dis)
        ind = 0
        # Loop over ocean_time steps:
        for n in range(len(ds.ocean_time.values)):
     
            if n == 0 or n == len(ds.ocean_time.values) -1:
                dt = 0.5* dtfac
            else:
                dt = dtfac
            print(n, ds.ocean_time.values[n], dt )

            ad_temp_n = np.zeros([ad_temp.shape[1], ad_temp.shape[2], ad_temp.shape[3]])
            ad_salt_n = np.zeros([ad_salt.shape[1], ad_salt.shape[2], ad_salt.shape[3]])


            temp_section = interpolator(Ginp, ds.isel(ocean_time = n), 'temp', lons, lats, ad_F=None)

            ad_ind_temporary = np.zeros_like(temp_section['val'])

            Km, Im = ad_ind_temporary.shape

            # Set all densities deeper than maxdepth  to zero
            mask = np.ones_like(temp_section['val'])
            mask[np.where( (temp_section['z_r'] - temp_section['dz']/2) < -1*maxdepth)] = 0

            for k in range(Km):
                    ad_ind_temporary[k, :] = ad_ind_temporary[k, :] + ad_ind

            # Set all values where mask = 0 to be 0 (simply multiply with mask as mask is either 1 or zero...). 
            ad_ind_temporary *= mask 
            ad_dens = -1* ad_ind_temporary * temp_section['z_r']*temp_section['dz']*temp_section['dislen']*1000 * dt

            ad_temp_array = -1* R0 * a * ad_dens
            ad_salt_array =  R0 * b * ad_dens

            # ADJOINT of interpolation here! 
            ad_temp_section = interpolator(Ginp, ds.isel(ocean_time = n), ad_temp_n, lons, lats, ad_F=ad_temp_array, mode = 'ad')
            ad_salt_section = interpolator(Ginp, ds.isel(ocean_time = n), ad_salt_n, lons, lats, ad_F=ad_salt_array, mode = 'ad')
            ad_temp[n, :, :, : ] = ad_temp_section['field']
            ad_salt[n, :, :, : ] = ad_salt_section['field']


        return ad_temp, ad_salt

def main():
    maxdepth = 98
    cyclenum = sys.argv[1]
    fileDir = sys.argv[2]
    constants_file = sys.argv[3]
    outputDir = sys.argv[4]

    force = 'false'
    if len(sys.argv) > 5:
        force = sys.argv[5]
    
    outputfile = outputDir + '/ocean_ads_outer0.nc-{}'.format(cyclenum).replace('//','/')
    fwdfile = fileDir + '/norshelf-ana_fwd_{0}T00Z_outer0.nc'.format(cyclenum).replace('//','/')
    #outputdir = /lustre/storeB/users/siljeci/FOCCUS
    #fileDir = /lustre/storeB/users/siljeci/FOCCUS
    #cyclenum = 20250501
    #fwdfile = fileDir + '/{0}/ocean_fwd_outer0.nc-{0}'.format(cyclenum).replace('//','/')

    if os.path.exists(outputfile) and force.lower() != 'true':
        print('File exists: {}'.format(outputfile))
        return
    
    ad_temp, ad_salt = PEA_forcing(cyclenum, fileDir, constants_file, maxdepth )

    adsfile = create_adsenfile(fwdfile, outputDir)
    
    with xr.open_dataset(fwdfile) as dsgrid:
        mask = dsgrid.mask_rho.values
    
    # Extend mask to same dimensions as ad_temp/ad_salt
    mask = np.broadcast_to(mask, (ad_temp.shape[0], ad_temp.shape[1],) + mask.shape)

    # Use 4D mask to set entries where mask == 0 to NaN
    ad_temp = np.where(mask == 0, np.nan, ad_temp)
    ad_salt = np.where(mask == 0, np.nan, ad_salt)

    ad_temp = np.ma.masked_invalid(ad_temp)
    ad_salt = np.ma.masked_invalid(ad_salt)

    with Dataset(adsfile, 'r+') as ds:
        np.ma.set_fill_value(ad_temp, ds.variables['v'][0,:].fill_value)
        np.ma.set_fill_value(ad_salt, ds.variables['v'][0,:].fill_value)

        ds.variables['temp'][:] = ad_temp
        ds.variables['salt'][:] = ad_salt

    print('Done with {}'.format(adsfile))
    command = 'ncks -7 -L 1 --ppc default=3#temp=8#salt=8  {} {}'.format(adsfile, adsfile.replace('_ads_', '_tmpads_'))
    subprocess.call(command, shell = True)
    command = 'mv {} {}'.format(adsfile.replace('_ads_', '_tmpads_'), adsfile)
    subprocess.call(command, shell = True)



if __name__ == "__main__":
    main()
