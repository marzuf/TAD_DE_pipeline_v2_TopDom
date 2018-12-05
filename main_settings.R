cat("load main_settings.R from TAD_DE_pipeline_v2_TopDom\n")

########################################################
### set paths to main files needed accross the workflow + general settings
########################################################

historyDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/data/gene_history_reordered.txt")

# file with mapping from entrez to chromosomic positions
#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.Rdata") # previous also held symbols
# holds only netrezID and positions
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")

# file with symbol synonyms to entrezID
#synoDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/all_entrez_syno_1.Rdata")
# mapping entrezID and all possible symbols
symbolDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_SYMBOL/final_entrez2syno.txt")

# mapping entrezID and ensemblID
ensemblDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_ENSEMBL/final_entrez2ensembl.txt")

# file with coordinates of all regions
TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    

# file with assignment from entrez to all regions
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

# match dataset with kind of process
process_file <- paste0(setDir, "/mnt/ed4/marie/other_datasets/dataset_process.csv")
process_file_with_samp <- paste0(setDir, "/mnt/ed4/marie/other_datasets/dataset_process_with_samp_annotated.csv")

nCpu <- 30

########################################################
### which kind of data to prepare and use across pipeline
########################################################

# if TRUE, from scripts 2_ to further, use only the genes that pass the DE analysis threshold for all analyses
useFilterCountData <- TRUE

# use only TADs (and corresponding genes)>= minSize et < maxQuantile and genes belonging to those TADs
useFilterSizeData <- TRUE

# if TRUE, from scripts 2_ to further, only focus on  TAD regions
useTADonly <- TRUE

########################################################
### Gene filter min CPM - for filter data (above) and 1_runGeneDE
########################################################

minCpmRatio <- 20/888 


########################################################
### TAD filter number of genes - for filter data (above) and used in 2_runWilcoxonTAD etc.
########################################################

# keep TAD with >= minNbrGeneTAD
minNbrGeneTAD <- 3

# keep TAD with <= maxQuant Genes !!! UPDATE HERE: FOR RUNNING PIPELINE V2 -> FOR TOPDOM, USE 0.99
maxQuantGeneTAD <- 0.99


########################################################
### permutation settings - 5_runPermutationsMedian.R        => STEP 5
########################################################

withExprClass <- TRUE

# number of permutations to run
nRandomPermut <- 10000

# number of classes of expression used for permutate expression data
permutExprClass <- 5

########################################################
### permutation settings (TAD permutation with fix size)- 5b_runPermutationsRandomTADsFixSize.R        => STEP 5b
########################################################
nRandomPermutFixSize <- 10
percentTADsizeFix <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5)
#percentTADsizeFix <- c(1)

########################################################
### permutation settings (TAD permutation with Gaussian dist. sizes)- 5c_runPermutationsRandomTADsGaussian.R        => STEP 5c
########################################################

nRandomPermutGaussian <- 100
percentTADsizeGaussian <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5)
#percentTADsizeGaussian <- c(1)

########################################################
### permutation settings (TAD permutation with shuffle)- 5d_runPermutationsRandomTADsShuffle.R        => STEP 5d
########################################################
nRandomPermutShuffle <- 10000

gene2tadAssignMethod <- "startPos"

########################################################
### ratio calculation also for the randomization of TADs => STEP 8
########################################################

step8_for_permutGenes <- TRUE
step8_for_randomTADsFix <- FALSE
step8_for_randomTADsGaussian <- FALSE
step8_for_randomTADsShuffle <- FALSE

########################################################
### number of TADs to select for plotting                   => STEP 13
########################################################

nTopLolliPlot <- 20
nTopVennPlot <- 50


########################################################
### for the plot with the area                              => STEP 14h
########################################################
step14_for_randomTADsFix <- FALSE
step14_for_randomTADsGaussian <- FALSE

########################################################
### for the plot with the area                              => STEP 14i
########################################################
step14_for_randomTADsShuffle <- FALSE

########################################################
### which ratios to compute                                 => STEPS 8c, 12c, 14c, 17c
########################################################

# allDown <- c("ratioDown", "FCdown", "meanRatioFC" , "prodConcord", "prodMeanConcord", "prodLogRatioNbr", "prodRatioSum", "prodSignedRatio")
# prodLogRatioNbr is not bounded => not used it anymore
# prodRatioSum is similar to prodSignedRatio
allDown <- c("ratioDown", "FCdown", "meanConcordRatioFC" , "prodConcord", "prodMeanConcord", "prodSignedRatio")
allDown <- c("ratioDown", "rescWeighted", "rescWeightedQQ", "prodSignedRatio")
#allDown <- "rescWeightedQQ"
allDown <- c("ratioDown",  "prodSignedRatio")

allDown_limited <- allDown

########################################################
### which quantile of the permutations to consider          => STEPS 8c, 12c, 14c, 17c
########################################################

permThresh <- 0.95

########################################################
### percent of the size of the random TADs compared to real TADs          => STEPS 8d
########################################################

percentTADsize <- 0.9

########################################################
### which scores to plant ranked                             => STEP 17
########################################################

# to build the table with % of TADs that are above a thresh 
#
ratioDown_thresh <- 0.2
prodSignedRatio_thresh <- 0.8

########################################################
### which scores to plant ranked                             => STEP 18
########################################################

# toPlotRanked <- c("ratioDown", "FCdown") in run settings NOW

########################################################
### signif threshold combined p-val                             => STEP 20
########################################################

adjCombinedPval_thresh <- 0.001

########################################################
### FDR threshold for empirical FDR                             => STEP 21
########################################################

empFDR_tresh  <- 0.20

########################################################
### thresholds functional enrichment analyses                             => STEP 23
########################################################

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


### settings for OUTPUT_FOLDER/*/23_funcEnrich_TopTableGenes_TADsGenes_v0:
nTopTADs <- 10
nTopTADs_Fisher <- 50
nTopGSEA <- 10
nTopGO <- 20
nTopKEGG <- 20
nPerm_gseGO <- 10000
nPerm_gseKEGG <- 10000
nPerm_GSEA <- 10000
p_adj_topGOthresh <- 0.05
p_adj_topKEGGthresh <- 0.05
p_adj_DE_thresh_Fisher <- 0.05

########################################################
### thresholds functional enrichment analyses                             => STEP 23topTADs
########################################################

### settings for OUTPUT_FOLDER/*/23topTADs_funcEnrich_TopTableGenes_TADsGenes_v0_20:
#funcTADs_nTopTADs_Fisher <- 20
#funcTADs_nTopTADs <- 20
#funcTADs_nPerm_gseGO <- 10000
#funcTADs_p_adj_DE_thresh_Fisher <- 0.05
#funcTADs_nShowPlot <- 10

### settings for OUTPUT_FOLDER/*/23topTADs_funcEnrich_TopTableGenes_TADsGenes:
funcTADs_nTopTADs_Fisher <- 50
funcTADs_nTopTADs <- 50
funcTADs_nPerm_gseGO <- 10000
funcTADs_p_adj_DE_thresh_Fisher <- 0.05
funcTADs_nShowPlot <- 10


########################################################
### thresholds functional enrichment analyses                             => STEP 23topTable
########################################################

funcTopTable_p_adj_DE_genes <- 0.05
funcTopTable_nPerm_gseGO <- 10000
funcTopTable_nShowPlot <- 10

########################################################
### thresholds semantic similarity                             => STEP 24
########################################################

semSimMetric <- "Wang"
combineSemSimMethod <- "BMA"
nTopTADs_semSim <- 10
nRandomSemSim <- 1000

########################################################
### to run without voom                             => STEP 1
########################################################

#runWithoutVoom <- c("GSE102073", "GSE71119")

