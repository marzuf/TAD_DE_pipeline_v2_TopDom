#!/usr/bin/bash

start_time=$(date -R)    

#set -e

settingFile="BUILDDT_SETTING_FILES/run_settings_buildDT.R"

# ./run_TAD_DE_TopDom_buildDTpval_v2.sh

#pipSteps=( "17c" "17e" "18" "18e" "18f" )
#pipSteps=( "17c" )
#pipSteps=( "18f" )

#pipSteps=( "18" )
#pipSteps=( "17f" )
#pipSteps=( "17fpre" "17hpre" "17ipre")
pipSteps=( "17ipreRatios")

outFold="BUILDDT_OUTPUT_FOLDER"


for step in "${pipSteps[@]}"; do
    echo "> START step $step"
    if [[ "$step" == "18f" ]] ; then

        echo Rscript ../TAD_DE_pipeline/18f_plotOrderedParam_limited.R $settingFile $outFold/18_plotOrderedParam/rank_scores_to_plot_DT.txt
        Rscript ../TAD_DE_pipeline/18f_plotOrderedParam_limited.R $settingFile $outFold/18_plotOrderedParam/rank_scores_to_plot_DT.txt

    else
        echo "./zzz_run_given_step_given_data_v2.sh $settingFile $step"
        ./zzz_run_given_step_given_data_v2.sh $settingFile $step
    fi
done




###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
