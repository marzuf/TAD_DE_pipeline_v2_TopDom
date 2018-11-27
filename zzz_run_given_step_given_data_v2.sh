#!/usr/bin/bash

# running like
#  ./zzz_run_given_step_given_data_v2.sh run_settings_TCGAcrc_entrez.R 14
#  ./zzz_run_given_step_given_data_v2.sh <run_file> <step#>

scriptFolder="../TAD_DE_pipeline_v2"

set -e

args=( "$@" )
args_len=${#args[@]}

settingF=${args[0]}

i=1
while [[ $i -lt args_len ]]; do
    scriptF="$(ls $scriptFolder | grep -P ^${args[$i]}_.+R$)"
    cmd="Rscript $scriptFolder/$scriptF $settingF"
    echo "> $cmd"
    $cmd
    ((i++))
done


