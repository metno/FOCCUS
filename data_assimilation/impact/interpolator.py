import roppy
import xarray as xr
import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset
import gsw
import toolbox

class SphereError(Exception):
    pass
class GridError(Exception):
    pass
class DimError(Exception):
    pass
class AdError(Exception):
    pass

def set_gridvars(gridT, ds):        
    if ds.spherical.values == 1:
        X = ds['lon_{}'.format(gridT)].values
        Y = ds['lat_{}'.format(gridT)].values
    else: 
        X = ds['x_{}'.format(gridT)].values
        Y = ds['y_{}'.format(gridT)].values

    if 'mask_{}'.format(gridT) in ds.variables.keys():
        mask = ds['mask_{}'.format(gridT)].values
    else: 
        mask = np.ones_like(X)
    return X , Y , mask


def isvector(v):
    if type(v) == list:
        v = np.array(v)

    shape = np.array(v.shape) 
    if shape.size == 1 and sum(shape) >1:
        return True
    elif shape.size == 1 and sum(shape) == 1:
        print('point')
        return False
    elif np.array(np.where(shape != 1)[0]).size > 1: 
        return False
    elif np.array(np.where(shape != 1)[0]).size == 1: 
        return True

def coordinate_matrix(grid):
    x   = np.linspace(0, grid.lat_rho.shape[1]-1, grid.lat_rho.shape[1])
    y   = np.linspace(0, grid.lat_rho.shape[0]-1, grid.lat_rho.shape[0])
    xi  = np.zeros_like(grid.lon_rho)
    yi  = np.zeros([grid.lon_rho.shape[1], grid.lon_rho.shape[0]])
    xi[:,:] = x
    yi[:,:] = y
    yi  = np.swapaxes(yi, 1, 0)
    return xi, yi

def interpolator(Ginp, ds, field, Xpath, Ypath, ad_F=None, mode ='fwd'):
    '''
    INTERPOLATOR: Spatially interpolate a ROMS field to horizontal
                  path/trajectory

    F = interpolator(Ginp, field, Xpath, Ypath, ad_F, [options])              

    Given any 2D or 3D ROMS variable FIELD, linearly interpolate to the
    requested path (Xpath, Ypath) in spherical or Cartesian coordinates.
    If input variable is 3D, the horizontally interpolated values are
    returned at all levels of the s-coordinate.
    
    Output includes elemental cell thicknesses and path segment lengths to
    facilitate spatial integration of the interpolated field.

    On Input: 

    Ginp   -  roppy Sgrid object
    ds     -  xarray dataset handle to ROMS output file (his/fwd/qck/avg).
    field  -  Name (string) of ROMS field to interpolate from (2D or 3D array). 
              Avoid pasing field time records. This function cannot be used
              for temporal interpolation.
    Xpath  -  Longitude (degrees_east) or Cartesian (m) coordinates in the 
              XI direction (1D or 2D array).
    Ypath  -  Latitude (degrees_north) or Cartesian (m) coordinates in the 
              ETA direction (1D or 2D array).
    [ad_F] -  Adjoint array of the same size of X if 2D array. If processing
              3D variables, it need an additional dimension for vertical 
              levels. Not used in the function, but added for compatability 
              with adjoint version of the interpolator function. 
    [mode] -  default is "normal forward interpolation". If  mode = 'ad' 
              the function will instead do the adjoint operation. 
              This requires that ad_F is given too. 
 
    % On Output:

    F           Interpolated field (dictionary)

                  F['field']      variable to process (2D or 3D)
                  F['Xgrd']       variable X-locations (degrees_east or m)
                  F['Ygrd']       variable Y-locations (degrees_north or m)
                  F['size']       variable size
                  F['Ctype']      variable C-grid type
                  F['method']     interpolation method
                  F['val']        interpolated values at (X,Y) or (X,Y,k)
                  F['Xval']       interpolation X-locations (input Xpath)
                  F['Yval']       interpolation Y-locations (input Ypath)
                  F['Lindex1']    linear indices (I  ,J  ) containing (X,Y)
                  F['Lindex2']    linear indices (I+1,J  ) containing (X,Y)
                  F['Lindex3']    linear indices (I+1,J+1) containing (X,Y)
                  F['Lindex4']    linear indices (I  ,J+1) containing (X,Y)
                  F['Jindex ']    lower-left field J-index containing (X,Y)
                  F['h']          bathymetry at requested (X,Y)
                  F['dis']        distance along the path (in km)
                  F['dislen']     path segment lengths for path integral
                  F['z_w']        cell center depths at requested (X,Y,N+1)
                  F['z_r']        cell face depths at requested (X,Y,N)
                  F['dz']         cell thicknesses at requested (X,Y,N)

    dis, z_r, dz and dislen are returned as 2D variables so that the
    interpolated output can be plotted with pcolor(F.dis,F.z_r,F.val)
    and integrated with sum(F.val .* F.dislen .* F.dz)

    Notice that this function is adjointed, see the function ad_interpolator.  It can
    be used to set-up functionals in various adjoint based algorithms, like
    adjoint sensitivities.

    The output linear indices are as follows:

                        Lindex4           Lindex3

                            4________________3
                            |                |
                            |                |
                            |<---p--->       |
                            |                |
                            |........o   :   |     field cell elements
                            |       .    :   |     containing (Xpath,Ypath)
                            |       .    q   |     as Matlab array
                            |       .    :   |
                            |____________:___|
                            1                2
                        Lindex1           Lindex2

            
    '''

    # Ensure constistency of input for the two different modes: 
    if not ad_F is None and mode == 'fwd':
        print('Warning: INTERPOLATOR does not use the ad_F input')
        print( 'Pass an empty argument')

    if mode == 'ad' and ad_F is None:
        raise AdError('AD_INTERPOLATOR: input argument ad_F is empty')
    


    # Initialize output structure.
    F = {'field'  : [], 'Xgrd'   : [], 'Ygrd'    : [],           
         'size'   : [], 'Cgrid'  : [], 'val'     : [], 
         'Lindex1': [], 'Lindex2': [], 'Lindex3' : [],
         'Lindex4': [], 'h'      : [], 
         'Xval'   : Xpath, 
         'Yval'   : Ypath 
         }
    
    # Set the interpolation method in scipy.interpolate griddata 
    F['method'] = 'linear'

    # Add fields related to trajectories/transects 
    if isvector(Xpath): 
        F['dis'] =  [] ; F['dislen'] = []
    else:
        xpath_y = Xpath.shape[0]
        xpath_x = Xpath.shape[1]

    # Remove possible singleton (time) dimension of input variable.
    if type(field) is str:
        F['field']  = np.squeeze(ds[field].values[:])
    else: 
        F['field']  =  np.squeeze(field)


    F['size'] = F['field'].shape

    if mode == 'ad' and not ad_F is None:
        # Check size of input ad_F
        Xprod    = np.prod(Xpath.shape)
        ad_Fprod = np.prod(ad_F.shape) 
        if len(F['field'].shape) == 3:
            N = F['field'].shape[0]
        else: 
            N = 1
        if ad_Fprod != Xprod*N: 
            print('\n**** Error AD_INTERPOLATOR - inconsistent size of array ad_F')
            print('\n\tprod(size(X))\t=  {}'.format(Xprod))
            print('\tprod(size(ad_F))\t= {}'.format(ad_Fprod))
            print('\tneeds to be \t = {}'.format(Xprod*N))
            print('\tsize(ad_F)\t= {}'.format(ad_F.shape))
            print('\tsize(field)\= {}'.format(field))
            raise AdError('Inconsistent size of array ad_F!!!')

        
    # Set ROMS grid size.
    # Lm  Number of I-direction INTERIOR RHO-points
    # Mm  Number of J-direction INTERIOR RHO-points
    L = F['field'].shape[-1] -1; Lp = F['field'].shape[-1]
    M = F['field'].shape[-2] -1; Mp = F['field'].shape[-2]

    is2d = False; is3d = False
    my_rank = len(F['field'].shape)

    if isvector(F['field']):
        raise DimError('Input FIELD is a 1D vector - that makes no sense')

    
    elif my_rank == 2:
        is2d = True
        Jm, Im = F['field'].shape

    elif my_rank == 3:
        is3d =  True 
        Km, Jm, Im = F['field'].shape
        # Add fields associated with 3D input variable to output structure
        F['z_w'] = []; F['z_r'] = []; F['dz'] = []; 
    
    else:
        raise DimError('Input FIELD must be a 2D or 3D array at 1 time slice.\nNumber of dimensions: {}'.format(my_rank))
    

    # Extract the ROMS grid lon/lat or x/y coordinates for this field.

    if (Im == L) & (Jm == M):
        F['Cgrid'] = 'psi_point'
        F['Xgrd'], F['Ygrd'], mask = set_gridvars(F['Cgrid'].split('_')[0], ds)    
        h = 0.25 * (Ginp.h[0:M, 0:L] + Ginp.h[1:Mp, 1:Lp] + 
                Ginp.h[0:M , 1:Lp] + Ginp.h[1:Mp, 0:L])
        
    elif (Im == Lp) & (Jm == Mp):
        F['Cgrid'] = 'rho_point'
        F['Xgrd'], F['Ygrd'], mask = set_gridvars(F['Cgrid'].split('_')[0], ds)    

        h = Ginp.h

    elif (Im == L) & (Jm == Mp):
        F['Cgrid'] = 'u_point'
        F['Xgrd'], F['Ygrd'], mask = set_gridvars(F['Cgrid'].split('_')[0], ds)    

        h = 0.5 * (Ginp.h[0:Mp, 0:L] + Ginp.h[0:Mp, 1:Lp] )

    elif (Im == Lp) & (Jm == M):
        F['Cgrid'] = 'v_point'
        F['Xgrd'], F['Ygrd'], mask = set_gridvars(F['Cgrid'].split('_')[0], ds)    

        h = 0.5 * (Ginp.h[0:M, 0:Lp] + Ginp.h[1:Mp, 0:Lp] )
    else:            
        raise GridError(' INTERPOLATOR: unable to determine Arakawa C-grid type:\n Cgrid is inconsistent with ROMS grid.')


    #--------------------------------------------------------------------------
    # Get bathymetry and z-coords along requested path.
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Set fractional coordinates and linear indices containing the point 
    # to interpolate. The fractional coordinates origin is irrelevant here
    # since is only used locally.
    #--------------------------------------------------------------------------
    xi, yi = coordinate_matrix(Ginp)

    xind = griddata( (F['Xgrd'].flatten(), F['Ygrd'].flatten()), xi.flatten(), (Xpath, Ypath) , method = F['method'] )
    yind = griddata( (F['Xgrd'].flatten(), F['Ygrd'].flatten()), yi.flatten(), (Xpath, Ypath) , method = F['method'] )
   
    xind1 = np.floor(xind).astype(int); xind2 = xind1 + 1; xind2[np.where(xind2 > Im)] = Im    # I, I+1
    yind1 = np.floor(yind).astype(int); yind2 =  yind1 + 1;  yind2[np.where(yind2 > Jm)] = Jm   # J, J+1

    # Interpweights: 
    p2 = (xind2 - xind1)*(xind - xind1); p1 = 1.-p2
    q2 = (yind2 - yind1)*(yind - yind1); q1 = 1.-q2

    w1 = mask[yind1, xind1] * p1*q1
    w2 = mask[yind1, xind2] * p2*q1
    w3 = mask[yind2, xind2] * p2*q2
    w4 = mask[yind2, xind1] * p1*q2
    mask_sum = mask[yind1, xind1] + mask[yind1, xind2] + mask[yind2, xind2] + mask[yind2, xind1]

    if any(mask_sum < 4):
        edit_inds = np.where(mask_sum < 4)[0]
        for edit_ind in edit_inds:
            wsum = w1[edit_ind] + w2[edit_ind] + w3[edit_ind] + w4[edit_ind]
            if wsum > 0:
                w1[edit_ind] /= wsum
                w2[edit_ind] /= wsum
                w3[edit_ind] /= wsum
                w4[edit_ind] /= wsum
            else: 
                w1[edit_ind] = 0
                w2[edit_ind] = 0
                w3[edit_ind] = 0
                w4[edit_ind] = 0

    F['Lindex1'] = [xind1, yind1]
    F['Lindex2'] = [xind2, yind1]
    F['Lindex3'] = [xind2, yind2]
    F['Lindex4'] = [xind1, yind2]

    F['h'] = w1*Ginp.h[yind1, xind1] + w2*Ginp.h[yind1, xind2] + w3*Ginp.h[yind2, xind2] + w4*Ginp.h[yind2, xind1] 

    if is3d:
        F['z_r'] = w1*Ginp.z_r[:,yind1, xind1] + w2*Ginp.z_r[:,yind1, xind2] + w3*Ginp.z_r[:,yind2, xind2] + w4*Ginp.z_r[:,yind2, xind1] 

        F['z_w'] = w1*Ginp.z_w[:,yind1, xind1] + w2*Ginp.z_w[:,yind1, xind2] + w3*Ginp.z_w[:,yind2, xind2] + w4*Ginp.z_w[:,yind2, xind1] 

        F['dz'] = np.diff(F['z_w'], n = 1, axis = 0 )
        
        if not isvector(Xpath):
            
            F['z_w'] = np.reshape(F['z_w'], newshape= (Km +1 , xpath_y, xpath_x))
            F['z_r'] = np.reshape(F['z_r'], newshape= (Km, xpath_y, xpath_x))
            F['dz'] = np.reshape(F['dz'], newshape= (Km, xpath_y, xpath_x))

    if not isvector(Xpath):
            F['h'] = np.reshape(F['h'], newshape= (xpath_y, xpath_x))
    


    #--------------------------------------------------------------------------
    # Path distance and segment lengths (for along-path integration).
    #--------------------------------------------------------------------------
    if isvector(Xpath):
        if ds.spherical.values == 1:
            # Use np.cumsum to calculate cumulative distance from s=0 is at the path points.
            # Use np.insert to ensure first entry = 0, and divide by 1000 to convert from 
            # meters to km.
            dis = np.insert(np.cumsum(gsw.distance(Xpath,Ypath)), 0, 0)/1000.

            # The cumulative distance from s=0 is at the path points. Average to
            # between-path points so that we can use diff to calculate path
            # segment lengths centered on the path points. But first Pad with the
            # first and last values of the original dis. Then diff will give half
            # cell widths for the end points of the path
            dtmp = 0.5 * (dis[:-1] + dis[1:])
            dtmp = np.append(np.insert(dtmp, 0, dis[0]),  dis[-1])
            dislen = np.diff(dtmp)
            if is2d:
                F['dis'] = dis
                F['dislen'] = dislen
            elif is3d:
                F['dis'] = np.repeat(np.reshape(dis, (1,len(Xpath))), Km, axis = 0)
                F['dislen'] = np.repeat(np.reshape(dislen, (1,len(Xpath))), Km, axis = 0)
        else:
            raise SphereError('Still need code for path in grid or normalized coordinates')
            
    else:
        print('Warning: No code yet to compute areas for input coordinates not a grid')

    #--------------------------------------------------------------------------
    # Interpolate. This is the only part that is adjointable.
    #--------------------------------------------------------------------------
    if is2d:
        # 2D field 
        F['val'] = w1*F['field'][F['Lindex1'][1], F['Lindex1'][0]] + w2*F['field'][F['Lindex2'][1], F['Lindex2'][0]] + w3*F['field'][F['Lindex3'][1], F['Lindex3'][0]] + w4*F['field'][F['Lindex4'][1], F['Lindex4'][0]] 

        if not isvector(Xpath):
                F['field'] = np.reshape(F['field'], newshape= (xpath_y, xpath_x))        
        if mode == 'ad': 
            # Adjoint of interpolation:
            for n in range(len(xind1)):
                F['field'][F['Lindex1'][1][n], F['Lindex1'][0][n]] = F['field'][F['Lindex1'][1][n], F['Lindex1'][0][n]] + w1[n]*ad_F[n]
                F['field'][F['Lindex2'][1][n], F['Lindex2'][0][n]] = F['field'][F['Lindex2'][1][n], F['Lindex2'][0][n]] + w2[n]*ad_F[n]
                F['field'][F['Lindex3'][1][n], F['Lindex3'][0][n]] = F['field'][F['Lindex3'][1][n], F['Lindex3'][0][n]] + w3[n]*ad_F[n]
                F['field'][F['Lindex4'][1][n], F['Lindex4'][0][n]] = F['field'][F['Lindex4'][1][n], F['Lindex4'][0][n]] + w4[n]*ad_F[n]
                ad_F[n] = 0
            
    elif is3d:
        F['val'] = w1*F['field'][:, F['Lindex1'][1], F['Lindex1'][0]] + w2*F['field'][:, F['Lindex2'][1], F['Lindex2'][0]] + w3*F['field'][:, F['Lindex3'][1], F['Lindex3'][0]] + w4*F['field'][:, F['Lindex4'][1], F['Lindex4'][0]] 

        if mode == 'ad':
            # Adjoint of interpolation:
            for n in range(len(xind1)):
                F['field'][:, F['Lindex1'][1][n], F['Lindex1'][0][n]] = F['field'][:, F['Lindex1'][1][n], F['Lindex1'][0][n]] + w1[n]*ad_F[:, n]
                F['field'][:, F['Lindex2'][1][n], F['Lindex2'][0][n]] = F['field'][:, F['Lindex2'][1][n], F['Lindex2'][0][n]] + w2[n]*ad_F[:, n]
                F['field'][:, F['Lindex3'][1][n], F['Lindex3'][0][n]] = F['field'][:, F['Lindex3'][1][n], F['Lindex3'][0][n]] + w3[n]*ad_F[:, n]
                F['field'][:, F['Lindex4'][1][n], F['Lindex4'][0][n]] = F['field'][:, F['Lindex4'][1][n], F['Lindex4'][0][n]] + w4[n]*ad_F[:, n]
                ad_F[:,n] = 0


        if not isvector(Xpath):
                F['field'] = np.reshape(F['field'], newshape= (Km, xpath_y, xpath_x))

    return F