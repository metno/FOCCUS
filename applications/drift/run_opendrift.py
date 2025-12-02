"""
Author: Mateusz Matuszak
"""

import sys
sys.path.insert(0, str("submodules/opendrift"))
from opendrift.models.oceandrift import OceanDrift

from datetime import datetime, timedelta
import os

def run_opendrift(file, lon, lat, z=0, N=1, radius=0, start_time=None, duration=12, time_step=30, time_step_output=60, outfile='sample_file.nc', depth_type='z', vertical_mixing=False, horizontal_diffusivity=0.1, coastline_shp=None):
    """
        A wrapper for running OpenDrift. https://opendrift.github.io/
    Args:
        file                    [str]       :   Model netCDF file containing ocean current data.    
        lon                     [float]     :   Initial longitude position.
        lat                     [float]     :   Initial latitude position.
        z                       [float]     :   Initial depth in meters.
        N                       [int]       :   Number of particles.
        radius                  [float]     :   Particle seeding radius around initial lon and lat. 
        start_time              [str]       :   Start time following %Y-%m-%dT%H:%M:%S format. 
        duration                [int]       :   Simulation duration in hours.
        time_step               [int]       :   Simulation time step in minutes. 
        time_step_output        [int]       :   Output frequency in minutes.
        outfile                 [str]       :   Name of output file. 
        depth_type              [str]       :   Vertical grid type of provided dataset. Currently supports 'z' and 's'. 
        vertical_mixing         [bool]      :   Set to True to allow particles to propagate vertically. False otherwise. 
        horizontal_diffusivity  [float]     :   Value for horizontal diffusivity. 
        coastline_shp           [str|list]  :   Add a custom coastline. Defaults to coastline from GSHHG.
    """
    #TODO add more tests for values
    if os.path.isdir(file):
        file = file+'/*'
    elif os.path.exists(file) or 'thredds' in file:
        pass
    else:
        raise ValueError(f'File or directory {file} not found.')
        
    if depth_type == 'z':
        from opendrift.readers import reader_netCDF_CF_generic
        r = reader_netCDF_CF_generic.Reader(file)
    elif depth_type == 's':
        from opendrift.readers import reader_ROMS_native
        r = reader_ROMS_native.Reader(file)
    else:
        raise ValueError(f'Supported depth types are [z, s], got {depth_type}')

    o = OceanDrift(
        loglevel=20,
        seed=0
    )
    o.add_reader(r)
    
    o.set_config('general:seafloor_action', 'lift_to_seafloor')
    o.set_config('general:coastline_action', 'previous')
    o.set_config('drift:horizontal_diffusivity', horizontal_diffusivity)
    o.set_config('drift:vertical_mixing', vertical_mixing)
    o.set_config('vertical_mixing:diffusivitymodel', 'environment')
    o.set_config('drift:advection_scheme', 'runge-kutta4')

    if coastline_shp is not None:
        if isinstance(coastline_shp, str):
            coastline_shp = [coastline_shp]
        if isinstance(coastline_shp, list) and all(isinstance(item, str) for item in coastline_shp):
            raise TypeError('Argument coastline_shp must be either a string or a list of strings.')
        else:
            from opendrift.readers import reader_shape
            for c in coastline_shp:
                r.append(reader_shape.Reader.from_shpfiles(c))

    if start_time is None:
        start_time = r.start_time
    elif type(start_time) == str:
        try:
            start_time = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
        except:
            raise ValueError('Provided start time does not follow required format. Try %Y-%m-%dT%H:%M:%S')
    elif type(start_time) == datetime.datetime:
        pass
    else:
        raise TypeError('Type of start_time is not supported.')
    
    o.seed_elements(lon=lon,
                    lat=lat,
                    z=z,
                    number=N,
                    radius=radius,
                    time=start_time)
    
    o.run(duration=timedelta(hours=duration),
          time_step=timedelta(minutes=time_step),
          time_step_output=timedelta(minutes=time_step_output),
          outfile=outfile,
          export_variables=[
              'sea_water_temperature',
              'sea_water_salinity',
              'z'
          ]
          )
    # for some reason it doesn't extract salinity and temp

if __name__ == "__main__":
    #TODO allow a list of lons, lats and z
    import argparse
    parser = argparse.ArgumentParser(
        prog='opendrift'
    )
    parser.add_argument(
        '-f', '--file', type=str, required=True, help='Input file or directory.'
    )
    parser.add_argument(
        '-lon', '--longitude', type=float, required=True, help='Initial longitude position of particles.'
    )
    parser.add_argument(
        '-lat', '--latitude', type=float, required=True, help='Initial latitude position of particles.'
    )
    parser.add_argument(
        '-z', '--depth_level', type=float, default=0, help='Intial depth level of particles in meters.'
    )
    parser.add_argument(
        '-N', '--number', type=int, default=1, help='Number of particles.'
    )
    parser.add_argument(
        '-r', '--radius', type=float, default=0, help='Radius for particle spread around initial lon lat position.'
    )
    parser.add_argument(
        '-st', '--start_time', type=str, default=None, help='Start time as ISO string. <year-month-dayThour:minutes:seconds.'
    )
    parser.add_argument(
        '-dur', '--duration', type=int, default=12, help='Duration of simulation in hours.'
    )
    parser.add_argument(
        '-ts', '--time_step', type=int, default=30, help='Computation time step during simulation in minutes.'
    )
    parser.add_argument(
        '-tso', '--time_step_output', type=int, default=60, help='Output frequency of simulation in minutes.'
    )
    parser.add_argument(
        '-o', '--outfile', type=str, default='sample_file.nc', help='Name of output file.'
    )
    parser.add_argument(
        '-dt', '--depth_type', default='z', choices=['z', 's'], type=str, help='Chose depth type of input data.'
    )
    parser.add_argument(
        '-vm', '--vertical_mixing', default=False, action=argparse.BooleanOptionalAction, help='When used, will turn on vertical mixing.'
    )
    parser.add_argument(
        '-hd', '--horizontal_diffusivity', default=0.1, type=float, help='Set horizontal diffusivity value.'
    )
    parser.add_argument(
        '-c' '--coastline_shp', default=None, help='Provide custom coastline from file(s).'
    )
    args = parser.parse_args()
    run_opendrift(file=args.file,
                  lon=args.longitude,
                  lat=args.latitude,
                  z=args.depth_level,
                  N=args.number,
                  radius=args.radius,
                  start_time=args.start_time, 
                  duration=args.duration,
                  time_step=args.time_step,
                  time_step_output=args.time_step_output,
                  outfile=args.outfile,
                  depth_type=args.depth_type,
                  vertical_mixing=args.vertical_mixing,
                  horizontal_diffusivity=args.horizontal_diffusivity,
                  coastline_shp=args.coastline_shp)

