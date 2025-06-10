#$ -S /bin/bash
#$ -l h_rt=02:00:00
#$ -q research-r8.q
#$ -l h_rss=48G,mem_free=48G,h_data=48G
#$ -o ./OUT_$JOB_NAME.$JOB_ID   #CHANGE THESE FOR LOGFILE DIRECTORY
#$ -e ./ERR_$JOB_NAME.$JOB_ID
#$ -N month-avr

# This script makes the MONTHLY average of YEAR

DIR_OUT= './' # SPECIFY OUTPUT DIR HERE

YEAR=2024

for MONTH in $(seq -w 1 12); do
    # Get the first and last day-of-year for the month
    first_day=$(date -d "2024-$MONTH-01" +%j)
    last_day=$(date -d "2024-$MONTH-01 +1 month -1 day" +%j)
    start_idx=$((10#$first_day - 1))
    end_idx=$((10#$last_day - 1))

    for OUTER in $(seq 0 1); do
        echo "Processing:" $YEAR $MONTH "outer:" $OUTER "first_day:" $first_day "last_day:" $last_day "start_idx:" $start_idx "end_idx:" $end_idx

        ncks -d time,$start_idx,$end_idx $DIR_OUT/norshelf-ana_fwd_avr_2024_outer${OUTER}.nc $DIR_OUT/tmp_${MONTH}_${OUTER}.nc
        echo "tmp_${MONTH}_${OUTER}.nc created, starting to average"
        ncwa -a time $DIR_OUT/tmp_${MONTH}_${OUTER}.nc $DIR_OUT/norshelf-ana_fwd_avr_2024${MONTH}_outer${OUTER}.nc
        #rm tmp_${MONTH}_${OUTER}.nc
    done
done

# TODO: something seems to be wrong with the time dim in the output files, in particular for december