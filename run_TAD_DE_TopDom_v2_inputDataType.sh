#!/usr/bin/bash

# ./run_TAD_DE_TopDom_v2_inputDataType.sh

start_time=$(date -R)    

#set -e

settingFold="../TAD_DE_pipeline/SETTING_FILES_cleanInput"

all_setting_files=($(ls $settingFold/*.R))


#all_setting_files=(
#"$settingFold/run_settings_GSE77314_normal_tumor.R"
#)

#all_setting_files=(
#"$settingFold/run_settings_GSE52166_prePf_postPf.R"
#)

#all_setting_files=(
#"$settingFold/run_settings_TCGApaad_paad_mutKRAS.R"
#"$settingFold/run_settings_TCGAacc_acc_mutCTNNB1.R"
#"$settingFold/run_settings_TCGAcrc_msi_mss.R"
#"$settingFold/run_settings_TCGAlaml_laml_mutFLT3.R"
#"$settingFold/run_settings_TCGAthca_thca_mutBRAF.R"
#"$settingFold/run_settings_TCGAlihc_lihc_mutCTNNB1.R"
#"$settingFold/run_settings_TCGAskcm_skcm_mutBRAF.R"
#"$settingFold/run_settings_TCGAluad_luad_mutKRAS.R"
#)

#all_setting_files=( "$settingFold/run_settings_TCGAskcm_skcm_mutCTNNB1.R" )

all_setting_files=( "$settingFold/run_settings_TCGAstad_EBVneg_EBVpos.R" )

#pipSteps=( "0cleanInput" "1cleanInput" )
# on electron:
# pipSteps=( "2" "3" "4" )
#pipSteps=( "5" "6" "7" "8c" )
#pipSteps=( "9" "10" "11" "13" "14f2" "14i2" "170v2" )
#pipSteps=( "19" "20" "21" "22" "23")
#pipSteps=( "8c" )
pipSteps=( "170" )

# on positron:
#pipSteps=( "5d" "14f2split" "14f2splitSlope" "14i2split" "14i2splitSlope" "24v2Wang" )

#pipSteps=( "170" )

#pipSteps=( "0cleanInput" "1cleanInput" "5" "5d" "8c" "14f2" "14i2")
#pipSteps=( "1cleanInput" )

#pipSteps=( "13cleanInput")
#pipSteps=( "5" "8c" "14f2" "14i2")
#pipSteps=( "14f2split" "14f2splitSlope" "14i2split" "14i2splitSlope" )
#pipSteps=( "19" )
#pipSteps=( "20" )
#pipSteps=( "21" )
#pipSteps=( "14f3" )
#pipSteps=( "14f4" )
#pipSteps=( "14i3" )
#pipSteps=( "14i4" )
#pipSteps=( "22" )
#pipSteps=( "170v2" )
#pipSteps=( "170prunePermGenes" )
#pipSteps=( "24v2" )
#pipSteps=( "23" )
#pipSteps=( "24vAllWang" )
#pipSteps=( "23topTADs" )
#pipSteps=( "23topTable" )

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
