#!/usr/bin/bash

start_time=$(date -R)    


#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES"

all_setting_files=($(ls $settingFold/*.R))


#pipSteps=( "0" "1" "2" "3" "4" "5" "5b" "5c" "5d" "6" "7" "8c" "9" "10" "11" "12c" "13" "14e" "14f" "14f2" "14h" "14i" )
#pipSteps=( "0" "1" "3" "4" "5" "5d" "8c" "14f2" "14i" )
#pipSteps=( "2" "6" "7" "9" "10" "11" "13" ) 
pipSteps=( "170" ) 

for settingF in "${all_setting_files[@]}"; do
    echo "> START for $settingF"
    for step in "${pipSteps[@]}"; do
        echo "./zzz_run_given_step_given_data_v2.sh $settingF $step"
        ./zzz_run_given_step_given_data_v2.sh $settingF $step
        #exit 0
    done 
    #exit 0
done




###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
