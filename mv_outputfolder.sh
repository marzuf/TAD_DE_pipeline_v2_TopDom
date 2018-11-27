#!/bin/bash

#for folder in OUTPUT_FOLDER/* ; do     
#    echo mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#    mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#done

#23_funcEnrich_TopTableGenes_TADsGenes

for folder in OUTPUT_FOLDER/*/23topTADs_funcEnrich_TopTableGenes_TADsGenes ; do     
    echo mv $folder ${folder}_v0_20
    mv $folder ${folder}_v0_20
#    echo mv $folder ${folder}_EXPRCLASS_NOT_FPKM
#    mv $folder ${folder}_EXPRCLASS_NOT_FPKM
done

#for folder in OUTPUT_FOLDER/*/23topTADs_funcEnrich_TopTableGenes_TADsGenes_v0_20 ; do     
#    echo mv $folder `dirname ${folder}`/23topTADs_funcEnrich_TopTableGenes_TADsGenes
#    mv $folder `dirname ${folder}`/23topTADs_funcEnrich_TopTableGenes_TADsGenes
#done


#for folder in OUTPUT_FOLDER/*/170_score_auc_pval_withShuffle ; do     
#    echo mv $folder `dirname ${folder}`/170v0_score_auc_pval_withShuffle
#    mv $folder `dirname ${folder}`/170v0_score_auc_pval_withShuffle
#done


### 23_funcEnrich_TopTableGenes_TADsGenes_v0
### settings for OUTPUT_FOLDER/*/23_funcEnrich_TopTableGenes_TADsGenes_v0:
#nTopTADs <- 10
#nTopTADs_Fisher <- 50
#nTopGSEA <- 10
#nTopGO <- 20
#nTopKEGG <- 20
#nPerm_gseGO <- 10000
#nPerm_gseKEGG <- 10000
#nPerm_GSEA <- 10000
#p_adj_topGOthresh <- 0.05
#p_adj_topKEGGthresh <- 0.05
#p_adj_DE_thresh_Fisher <- 0.05

### 23_funcEnrich_TopTableGenes_TADsGenes_v1
### settings for OUTPUT_FOLDER/*/23_funcEnrich_TopTableGenes_TADsGenes_v1:
#nTopTADs <- 20
#nTopTADs_Fisher <- 100
#nTopGSEA <- 20
#nTopGO <- 50
#nTopKEGG <- 50
#nPerm_gseGO <- 10000
#nPerm_gseKEGG <- 10000
#nPerm_GSEA <- 10000
#p_adj_topGOthresh <- 0.05
#p_adj_topKEGGthresh <- 0.05
#p_adj_DE_thresh_Fisher <- 0.05


### 23_funcEnrich_TopTableGenes_TADsGenes_v2
### settings for OUTPUT_FOLDER/*/23_funcEnrich_TopTableGenes_TADsGenes_v2:
#nTopTADs <- 50
#nTopTADs_Fisher <- 150
#nTopGSEA <- 50
#nTopGO <- 100
#nTopKEGG <- 100
#nPerm_gseGO <- 10000
#nPerm_gseKEGG <- 10000
#nPerm_GSEA <- 10000
#p_adj_topGOthresh <- 0.05
#p_adj_topKEGGthresh <- 0.05
#p_adj_DE_thresh_Fisher <- 0.05
