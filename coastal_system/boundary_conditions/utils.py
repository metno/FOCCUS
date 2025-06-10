import matplotlib.pyplot as plt
import numpy as np

figures = []
figure_counter = 1

def plot_2boundary_slice(var_name, var_west, var_west_obc, **kwargs):
    global figure_counter
    fig, axs = plt.subplots(1, 2, sharey=True, sharex=True)
    (ax1, ax2) = axs
    fig.suptitle(f'West boundary - {var_name}')
    ax1.set_xlabel("xi_rho")
    ax2.set_xlabel("xi_rho")
    ax1.set_ylabel("s_rho")
    fig.subplots_adjust(wspace=0.03, left=0.06, right=1, bottom=0.1)
    #fig.subplots_adjust(wspace=0.03)

    im1 = ax1.pcolor(var_west,**kwargs)
    im2 = ax2.pcolor(var_west_obc,**kwargs)
    
    fig.colorbar(im2, ax=axs, pad=0.01, label=f'{var_name}')

    figures.append((f'figure_{figure_counter}', fig))
    figure_counter += 1

    #TODO:
    # should also plot countur lines in figure:
    # , cvar_west=None, cvar_obc_west=None, clevels=None
    #cn1 = ax1.contour(cvar_west, linestyles='--', cmap='Reds', levels=clevels)
    #cn2 = ax2.contour(cvar_obc_west, linestyles='--', cmap='Reds', levels=clevels) 

    return fig, axs

def get_boundary_slice(ds,var_name,direction=None,itime=0,obc=False,obc_adjust=0):
    # Returns an array of length 901 no matter which boundary
    # If itime is an array several time steps are returned. Default to return only first timestep
    
    if direction==None or direction not in ['west', 'south', 'east', 'north']:
        raise KeyError("You should give the direction of the boundary as input! (west,south,east,north)")

    if obc:
        # Variables with ending "_obc" contains 2D boundaries / slices (+ time dimension)
        # which are selected from the index 'boundary'
        ib = {'south': 0, 'east': 1, 'north': 2, 'west': 3}
        var = ds[var_name+'_obc'].isel(ocean_time=itime, boundary=ib[direction], obc_adjust=obc_adjust).values

    else:
        # Select the boundary slice from the 4D data
        # First get the variable:
        var_itime = ds[var_name].isel(ocean_time=itime).values

        # The indices in x and y-direction depending on the variable
        if var_name in ['u','ubar','u_obc','ubar_obc']:
            n_x = 900
            n_y = 351
        elif var_name in ['v','vbar','v_obc','vbar_obc']:
            n_x = 901
            n_y = 350
        else: 
            n_x = 901
            n_y = 351

        if isinstance(itime, int):
            # Get one time step
            var = np.full((42, 901), np.nan)
            if direction == 'south':
                var[:, :n_y] = var_itime[:, :, 0] 
            elif direction == 'east':
                var[:, :n_x] = var_itime[:, 0, :]
            elif direction == 'north':
                var[:, :n_y] = var_itime[:, :, -1]
            elif direction == 'west':
                var[:, :n_x] = var_itime[:, -1, :]
        elif isinstance(itime, np.ndarray):
            # Get many time steps
            var = np.full((len(itime), 42, 901), np.nan)
            if direction == 'south':
                var[:, :, :n_y] = var_itime[:, :, :, 0]
            elif direction == 'east':
                var[:, :, :n_x] = var_itime[:, :, 0, :]
            elif direction == 'north':
                var[:, :, :n_y] = var_itime[:, :, :, -1]
            elif direction == 'west':
                var[:, :, :n_x] = var_itime[:, :, -1, :]

    return var


def plot_domain_maplayer(ax,var):
    #TODO
    ax.pcolor(var) # or pcolormesh


# If want to scatter-plot need this:
#x = ds0.xi_rho.values
#y = ds0.eta_rho.values
#xv, yv = np.meshgrid(x, y)