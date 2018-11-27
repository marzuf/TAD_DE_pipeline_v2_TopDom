#!/usr/bin/bash

# running like
#  ./zzz_run_given_step_given_data.sh run_settings_TCGAcrc_entrez.R 14
#  ./zzz_run_given_step_given_data.sh <run_file> <step#>

maxJobs=20
maxLoad=60
Rexec="Rscript"
scriptFolder="../TAD_DE_pipeline_v2"

# ./zzz_run_given_step_given_data.sh run_settings_GSE40419_normal_cancer.R 0

settingFolder="../TAD_DE_pipeline/SETTING_FILES_cleanInput"

set -e

#step="23vAll"
#step="23"
#step="23topTable"
step="23topTADs"
#all_settingF=(
#"run_settings_GSE48166_control_ICM.R"
##"run_settings_GSE73765_noninf_list.R"
##"run_settings_GSE81046_salm_list.R"
#)


all_settingF=( $(ls $settingFolder/*.R))

#for sF in ${all_settingF[@]}; do
#echo $sF
#done
#exit

#i=1
#while [[ $i -lt args_len ]]; do
#    scriptF="$(ls $scriptFolder | grep -P ^${step}_.+R$)"
#    cmd="Rscript $scriptFolder/$scriptF $settingF"
#    echo "> $cmd"
#    $cmd
#    ((i++))
#done


#   scriptF="$(ls $scriptFolder | grep -P ^${step}_.+R$)"; cmd="Rscript $scriptFolder/$scriptF $settingF"; echo "> $cmd"; $cmd

#parallel -i -j $maxJobs -l $maxLoad sh -c "scriptF='HELLO'; echo $scriptF" -- ${all_settingF[@]}

#   scriptF="$(ls $scriptFolder | grep -P ^${step}_.+R$)"; cmd="Rscript $scriptFolder/$scriptF $settingF"; echo "> $cmd"; $cmd

#parallel -i -j $maxJobs -l $maxLoad sh -c 'echo '$Rexec' '$scriptFolder'/$(ls '$scriptFolder' | grep -P ^'${step}'_.+R$) '$settingFolder'/'{}' ; 
#                                           '$Rexec' '$scriptFolder'/$(ls '$scriptFolder' | grep -P ^'${step}'_.+R$) '$settingFolder'/'{} -- ${all_settingF[@]}

parallel -i -j $maxJobs -l $maxLoad sh -c 'echo '$Rexec' '$scriptFolder'/$(ls '$scriptFolder' | grep -P ^'${step}'_.+R$) '{}' ; 
                                           '$Rexec' '$scriptFolder'/$(ls '$scriptFolder' | grep -P ^'${step}'_.+R$) '{} -- ${all_settingF[@]}



#parallel -i -j $maxJobs -l $maxLoad sh -c 'cmd='$Rexec'; echo '$cmd-- ${all_settingF[@]}

# ./zzz_run_given_step_given_data.sh BUILDDT_SETTING_FILES/run_settings_buildDT.R

#parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript compare_chromo.R {}" -- ${all_settingF[@]}

#    parallel -i -j $maxJobs -l $maxLoad sh -c $Rexec' '$step1_script' -f $(ls '$step0outfolder'/'${cell_line}'_'${treatment}'_'{}'_'${binSize}'_aggregCounts.txt) -c '{}' -n '$norm_meth' -o '$step1outfolder' -b '$binSize' -t '$norm_caller -- ${chromo[@]}










