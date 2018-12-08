#!/usr/bin/bash

start_time=$(date -R)    

#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES_cleanInput"

#all_setting_files=($(ls $settingFold/*.R))

# had 0 column sample
#GSE52166_prePf_postPf.R

# error DEseq2
#GSE64810_control_carrier.R
#GSE101521_control_mdd.R
#GSE102073_stic_nostic.R
all_setting_files=(
"$settingFold/run_settings_GSE52166_prePf_postPf.R"
"$settingFold/run_settings_GSE64810_control_carrier.R"
"$settingFold/run_settings_GSE101521_control_mdd.R"
"$settingFold/run_settings_GSE102073_stic_nostic.R"
)


#pipSteps=( "0" "1" "2" "3" "4" "5" "5b" "5c" "5d" "6" "7" "8c" "9" "10" "11" "12c" "13" "14e" "14f" "14f2" "14h" "14i" )
#pipSteps=( "5b" "5c" "8c" "14h" )
#pipSteps=( "5d" "8c" "14h" "14i" )
pipSteps=( "0cleanInput" "1cleanInput" "5" "5d" "8c" "14f2" "14i2")
#pipSteps=( "1cleanInput" )

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
