#!/usr/bin/bash

start_time=$(date -R)    


all_setting_files=(
"OUTPUT_FOLDER/GSE101521_control_mdd"
"OUTPUT_FOLDER/GSE102073_stic_nostic"
"OUTPUT_FOLDER/GSE48166_control_ICM"
"OUTPUT_FOLDER/GSE51799_control_carrier"
"OUTPUT_FOLDER/GSE52166_prePf_postPf"
"OUTPUT_FOLDER/GSE57148_normal_COPD"
"OUTPUT_FOLDER/GSE58135_ERpos_adjERpos"
"OUTPUT_FOLDER/GSE58135_ERpos_tripleNeg"
"OUTPUT_FOLDER/GSE58135_tripleNeg_adjTripleNeg"
"OUTPUT_FOLDER/GSE61476_unaffectNPCday11_unaffectNPCday31"
"OUTPUT_FOLDER/GSE64810_control_carrier"
"OUTPUT_FOLDER/GSE64813_prePTSD_postPTSD"
"OUTPUT_FOLDER/GSE65540_before_after"
"OUTPUT_FOLDER/GSE66306_before_after"
"OUTPUT_FOLDER/GSE67528_ASDiPSC_ASDneuron"
"OUTPUT_FOLDER/GSE67528_ASDiPSC_ASDnpc"
"OUTPUT_FOLDER/GSE67528_CTRLiPSC_CTRLneuron"
"OUTPUT_FOLDER/GSE67528_CTRLiPSC_CTRLnpc"
)


for currfile in "${all_setting_files[@]}"; do
    ls $currfile | wc -l
done

#exit

#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES"

all_setting_files=(
"$settingFold/run_settings_GSE101521_control_mdd.R"
"$settingFold/run_settings_GSE102073_stic_nostic.R"
"$settingFold/run_settings_GSE48166_control_ICM.R"
"$settingFold/run_settings_GSE51799_control_carrier.R"
"$settingFold/run_settings_GSE52166_prePf_postPf.R"
"$settingFold/run_settings_GSE57148_normal_COPD.R"
"$settingFold/run_settings_GSE58135_ERpos_adjERpos.R"
"$settingFold/run_settings_GSE58135_ERpos_tripleNeg.R"
"$settingFold/run_settings_GSE58135_tripleNeg_adjTripleNeg.R"
"$settingFold/run_settings_GSE61476_unaffectNPCday11_unaffectNPCday31.R"
"$settingFold/run_settings_GSE64810_control_carrier.R"
"$settingFold/run_settings_GSE64813_prePTSD_postPTSD.R"
"$settingFold/run_settings_GSE65540_before_after.R"
"$settingFold/run_settings_GSE66306_before_after.R"
"$settingFold/run_settings_GSE67528_ASDiPSC_ASDneuron.R"
"$settingFold/run_settings_GSE67528_ASDiPSC_ASDnpc.R"
"$settingFold/run_settings_GSE67528_CTRLiPSC_CTRLneuron.R"
"$settingFold/run_settings_GSE67528_CTRLiPSC_CTRLnpc.R"
)




#pipSteps=( "0" "1" "2" "3" "4" "5" "5b" "5c" "5d" "6" "7" "8c" "9" "10" "11" "12c" "13" "14e" "14f" "14f2" "14h" "14i" )
pipSteps=( "5d" "8c" "14i" )



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
