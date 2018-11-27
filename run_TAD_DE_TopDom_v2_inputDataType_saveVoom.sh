#!/usr/bin/bash

start_time=$(date -R)    

#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES_cleanInput"

all_setting_files=($(ls $settingFold/*.R))

#all_setting_files=(
#"$settingFold/run_settings_TCGAcrc_msi_mss.R"
#"$settingFold/run_settings_GSE71119_dediffSM_MFSM.R"
#)

pipSteps=( "1saveVoomOnly" )

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
