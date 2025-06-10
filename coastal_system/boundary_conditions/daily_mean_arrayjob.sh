#$ -S /bin/bash
#$ -l h_rt=2:00:00
#$ -q research-r8.q
#$ -l h_rss=24G,mem_free=24G,h_data=24G
#$ -o ./OUT_$JOB_NAME.$JOB_ID   #CHANGE THESE FOR LOGFILE DIRECTORY
#$ -e ./ERR_$JOB_NAME.$JOB_ID
#$ -N day-avr-array
#$ -t 1-366

# Array jobs: This script makes the daily average for X dates after start_date. One job per date/file.

DIR_NORSHELF='https://thredds.met.no/thredds/dodsC/foccus/'
DIR_OUT= './' # SPECIFY OUTPUT DIR HERE

start_date=20240101

# Function to increment date
increment_date() {
    date -d "$1 + 1 day" +"%Y%m%d"
}

# Get a list of all dates to select the correct date for the array job
date_now=$start_date
all_dates=()
all_dates+=("placeholder_first_element") # since $SGE_TASK_ID starts counting from 1
for i in $(seq 1 365);
do
    all_dates+=("$date_now")
    date_now=$(increment_date "$date_now")
done

# get the date for this job
job_date="${all_dates[$SGE_TASK_ID]}"

YEAR=${job_date:0:4}
MONTH=${job_date:4:2}
DAY=${job_date:6:2}

for OUTER in $(seq 0 1);
do
    echo "Processing:" $job_date "outer:" $OUTER #$YEAR $MONTH $DAY
    ncwa -a ocean_time $DIR_NORSHELF/$YEAR/$MONTH/norshelf-ana_fwd_${YEAR}${MONTH}${DAY}T00Z_outer${OUTER}.nc $DIR_OUT/norshelf-ana_fwd_avr_${YEAR}${MONTH}${DAY}_outer${OUTER}.nc
done

# then do 
# ./concat_add_time_dim.sh

 