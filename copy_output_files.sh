#!/bin/bash

mainFold="OUTPUT_FOLDER"

#outFold="tmp_svg_14f2"
#mkdir -p $outFold

#all_folders=( $( ls $mainFold ) )

#for folder in "${all_folders[@]}" ; do     
#    echo $folder
#    cp $mainFold/$folder/14f2_cumulAllDown_limited_AUC/allRatios_cumsum_obs_permut.svg $outFold/${folder}_allratios_permGenes.svg
#    #echo mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#    # mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#done



#outFold="tmp_svg_19"
#mkdir -p $outFold

#all_folders=( $( ls $mainFold ) )

#for folder in "${all_folders[@]}" ; do     
#    echo $folder
#    cp $mainFold/$folder/19_SAM_emp_measurement/FDR_var_cut_off_logFC.svg $outFold/${folder}_FDR_var_cut_off_logFC.svg
#    cp $mainFold/$folder/19_SAM_emp_measurement/FDR_var_cut_off_intraTADcorr.svg $outFold/${folder}_FDR_var_cut_off_intraTADcorr.svg
#    cp $mainFold/$folder/19_SAM_emp_measurement/FDR_var_cut_off_ratioDown.svg $outFold/${folder}_FDR_var_cut_off_ratioDown.svg
#    cp $mainFold/$folder/19_SAM_emp_measurement/FDR_var_cut_off_prodSignedRatio.svg $outFold/${folder}_FDR_var_cut_off_prodSignedRatio.svg
#    #echo mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#    # mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#done


#outFold="tmp_svg_14f3"
#mkdir -p $outFold

#all_folders=( $( ls $mainFold ) )

#for folder in "${all_folders[@]}" ; do     
#    echo $folder
#    cp $mainFold/$folder/14f3_cumulAllDown_withDist_AUC/allRatios_cumsum_obs_permut_withDist.svg $outFold/${folder}_allRatios_cumsum_obs_permut_withDist.svg
##    mv $mainFold/$folder/14f3_cumulAllDown_withDist_AUC/allRatios_cumsum_obs_permut.svg $mainFold/$folder/14f3_cumulAllDown_withDist_AUC/allRatios_cumsum_obs_permut_withDist.svg
#done

all_folders=( $( ls $mainFold ) )

for folder in "${all_folders[@]}" ; do     
    echo $folder
    mv $mainFold/$folder/170prunePermGenes_score_auc_pval_withShuffle $mainFold/${folder}/170prunePermGenes_score_auc_pval_withShuffle_i50

done





