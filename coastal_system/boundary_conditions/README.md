This directory is intended for code relating to improvements in the model interfaces between CMEMS and MSCS. 

# How to use this directory

The base data files are (ex. for 2024-03-10):

`norshelf-ana_fwd_20240310T00Z_outer*.nc` 

where `*` should be replaced with `0` or `1`. These files contains 17 time steps and many variables. We are interested in the variables with a name ending `_obc`, which only have values at the boundary, and therefore have a different grid from other "normal" variables.  See the files `timestep_selection.ipynb` and `data_explotarion.ipynb` for more details about the variables, their dimensions and the how the time steps in the files overlap. 

## Initial data analysis

In order to reduce the amount of data the initial analysis, a monthly mean of the above files have been obtained for 2024. Below is some information about the post-processing performed to obtain first the daily, then the montly means. 

### Post-processing / obtaining mean files

1. Run `daily_mean_arrayjob.sh` to use `ncwa` to average over all hours of each day in a year. That is, this file loops over all the files in 2024 and their 17 time steps and outputs an average of all time steps for each day giving 366 files for 2024.

2. Run `concat_add_time_dim.sh`. This script takes the .nc-files averaged in the above step and adds a time dimension, sets the unit for this dimension and creates a record dimension and writes this to temporary files. Then these temorary files are concatenated into one file for 2024, which then contains 366 time steps. 

3. Run `monthly_mean.sh`. This script does step 1. and 2. above, but chreates monthly means. It takes the daily average files outputted in step 2. and creates a file with 12 time steps. 

4. Run `get_monthly_zdepth.py`. This file takes the files from step 3. and uses the `xroms` package to transform the data so that we may easily plot the depth in meters. The outputs two files: `norshelf-ana_fwd_avr_zdepth_202401-12_outer*.nc` where `*` is`0` or `1`.

## Plotting

Main script: `bias.ipynb`
**TODO: change path to thredds**

See the difference of `salt_obc_outer0 - salt_obc_outer1` or `temp_obc_outer0 - temp_obc_outer1` per month, or yearly for the western boundary. 

This script may easily be expanded to plot the other boundaries or other variables (as long as they were included in the above post-processing).



## TODO - Add more details below

Only files to show on GitHub will be `bias.ipynb` and the threadds file???

### DATAFILES

- change output dir in the above scripts --> foccus?
- create symlink for `norshelf-ana_fwd_avr_zdepth_202401-12_outer*.nc` to aviod double storage

### Information files

Use these files to understand the structure of the datafiles. 

- `timestep_selection.ipynb`
- `data_explotarion.ipynb`
- `compare_bry.ipynb` (unfinished, should use z-depth when comparing?)

### Utility files

- `utils.py`: contains scripts for selecting boundary data for `obc` variables or all.
- `check_diff.ipynb`: Quickly see the magnitude of `obc`-diff. Simply output the difference between `obc` outer 0 and 1 for a specific date for each month. 



