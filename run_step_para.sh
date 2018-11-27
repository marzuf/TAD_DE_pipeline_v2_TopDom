#!/usr/bin/bash

start_time=$(date -R)    

#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES_cleanInput"

all_setting_files=($(ls $settingFold/*.R))

maxJobs=30
maxLoad=40


parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript ../TAD_DE_pipeline_v2/19shuffle_SAM_emp_measurement.R {}" -- ${all_setting_files[@]}

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
