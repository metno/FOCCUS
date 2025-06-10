#$ -S /bin/bash
#$ -l h_rt=01:00:00
#$ -q research-r8.q
#$ -l h_rss=4G,mem_free=4G,h_data=4G
#$ -o ./OUT_$JOB_NAME.$JOB_ID   #CHANGE THESE FOR LOGFILE DIRECTORY
#$ -e ./ERR_$JOB_NAME.$JOB_ID
#$ -N add-time-dim-concat

# This script takes the output of 
# daily_mean_arrayjob.sh
# and adds the correct netcdf time dims. etc., then concats the induvidual 
# files to one larger netcdf file for a whole year. 
# NOTE: remove all intermediate files manually upon successful completion. 

DIR_OUT= './' # SPECIFY OUTPUT DIR HERE

start_date=20240101

# Function to increment date 
increment_date() {
    date -d "$1 + 1 day" +"%Y%m%d"
}

# Get a list of all dates to select the correct date (handle missing dates)
date_now=$start_date
for i in $(seq 1 366);
do
    YEAR=${date_now:0:4}
    MONTH=${date_now:4:2}
    DAY=${date_now:6:2}
    i_day=$((i-1))

    for OUTER in $(seq 0 1);
    do
        file=$DIR_OUT/norshelf-ana_fwd_avr_${YEAR}${MONTH}${DAY}_outer${OUTER}.nc
        if [ -f $file ]; then
            outfile=${file}_tmp
            echo "Processing:" $date_now "outer:" $OUTER $i #$YEAR $MONTH $DAY $file

            # set new dim and value
            ncap2 -O -s "time=${i_day}" $file $outfile
            # set unit of dim
            ncatted -O -a units,time,c,c,"days since ${YEAR}-01-01" $outfile
            # create record dim
            ncecat -O -u time $outfile $outfile
        fi
    done
    # get the next date
    date_now=$(increment_date "$date_now")
done
echo "Induvidual file processing done, concatenating *.nc_tmp files" $YEAR 

# Concatenate the files
for OUTER in $(seq 0 1);
    do
    ncrcat -O -d time,0,366 $DIR_OUT/*outer${OUTER}.nc_tmp $DIR_OUT/norshelf-ana_fwd_avr_${YEAR}_outer${OUTER}.nc
done

# check results:
ncdump -v time $DIR_OUT/norshelf-ana_fwd_avr_${YEAR}_outer0.nc
ncdump -v time $DIR_OUT/norshelf-ana_fwd_avr_${YEAR}_outer1.nc

# After finishing you may manually remove all intermediate files 
# and keep only $DIR_OUT/norshelf-ana_fwd_avr_${YEAR}_outer${OUTER}.nc 
