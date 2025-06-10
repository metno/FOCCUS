"""
Example runs of OpenDrift trajectory simulations
Author: Mateusz Matuszak
"""

from run_opendrift import run_opendrift

run_opendrift(
    file = 'https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_160m_m70_be',
    lon = 6.19,
    lat = 62.384,
    start_time='2024-06-12T12:00:00',
    duration=250,
    time_step=60,
    outfile='sample_file.nc',
    horizontal_diffusivity=0
)
