startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GOsignifEntropy.R\n"))

# Rscript cmp_DE_TADs_genes_GOsignifEntropy.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

buildTable <- TRUE

library(reshape2)
library(ggplot2)

source("analysis_utils.R")

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

registerDoMC(ifelse(SSHFS, 2, 30))

topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

outFold <- "CMP_TADs_GENES_GO_SIGNIF_ENTROPY"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_signif_logFile.txt")
system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

# GO for BP nad MF [do not take c5_CC]
gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/c5.bp_mf.v6.1.entrez.gmt")
stopifnot(file.exists(gmtFile))

# settings for enrichment analysis: clusterProfiler::enricher
# default values enricher function: (universe as background !)
# pvalueCutoff = 0.05 
# pAdjustMethod = "BH"
# minGSSize = 10
# maxGSSize = 500
# qvalueCutoff = 0.2 
enricher_pvalueCutoff <- 1
enricher_pAdjustMethod <- "BH"
enricher_minGSSize <- 1
enricher_maxGSSize <- 500
enricher_qvalueCutoff <- 1

enricher_results_sortGOby <- "p.adjust"
enricher_results_cmpNbrTop <- 50

mySub <- paste0("(# top GO: ",enricher_results_cmpNbrTop, ")") 

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... gmtFile:\t", gmtFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pvalueCutoff:\t", enricher_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pAdjustMethod:\t", enricher_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher minGSSize:\t", enricher_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher maxGSSize:\t", enricher_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher qvalueCutoff:\t", enricher_qvalueCutoff, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:10]

c5_msigdb <- read.gmt(gmtFile)

# all_ds <- all_ds[1:3]
curr_ds = all_ds[1]
if(buildTable){
  all_ds_DT <- foreach(curr_ds = all_ds, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)
    
    ### RETRIEVE TAD PVALUES
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    tad_pval_sort <- sort(tad_pval, decreasing=F)
    tad_pval_rank <- rank(tad_pval, ties.method="min")
    tad_pval_rank <- sort(tad_pval_rank)
    stopifnot( names(tad_pval_rank[tad_pval_rank==1])[1] == names(tad_pval_sort[1]) ) 
    
    nTopTADs <- round(topPercent * length(tad_pval_sort))
    stopifnot(nTopTADs > 0)
    
    topTADs <- names(tad_pval_sort[1:nTopTADs])
    
    ### RETRIEVE GENES PVALUES
    step1_fold <- file.path(dsFold, curr_ds, "1_runGeneDE")
    stopifnot(file.exists(step1_fold))
    toptableFile <- file.path(step1_fold, "DE_topTable.Rdata")
    stopifnot(file.exists(toptableFile))
    topTable_DT <- eval(parse(text = load(toptableFile)))
    stopifnot(!any(duplicated(topTable_DT$genes)))
    stopifnot(topTable_DT$genes == rownames(topTable_DT))
    gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(names(pipelineGenes) %in% topTable_DT$genes)
    
    topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipelineGenes),]
    genes_entrez <- sapply(topTable_DT$genes, function(x) pipelineGenes[x])
    
    stopifnot(!is.na(genes_entrez))
    stopifnot(length(genes_entrez) == nrow(topTable_DT) )
    stopifnot(names(pipelineGenes) %in% names(gene_pval))
    # gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    entrez_pval <- setNames(topTable_DT$adj.P.Val, genes_entrez)
    entrez_pval_sort <- sort(entrez_pval, decreasing=F)
    entrez_pval_rank <- rank(entrez_pval, ties.method="min")
    
    stopifnot( names(entrez_pval_rank[entrez_pval_rank==1])[1] == names(entrez_pval_sort[1]) ) 
    stopifnot( names(entrez_pval_rank) %in% gene2tad_DT$entrezID )
    stopifnot( names(tad_pval_rank) %in% gene2tad_DT$region )
    
    nTopGenes <- round(topPercent * length(entrez_pval_sort))
    # take the top percent DE genes
    topGenes_manyAsTopPercent <- names(entrez_pval_sort[1:nTopGenes])
    
    stopifnot(nTopGenes > 0)
    stopifnot(topTADs %in% gene2tad_DT$region)
    
    topTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% topTADs]  
    topTADs_genes <- topTADs_genes[topTADs_genes %in% pipelineGenes]
    stopifnot(length(topTADs_genes) > 0)
    stopifnot(length(topTADs_genes) <= length(entrez_pval_sort) )
  
    # take as many top genes as genes from topTADs  
    topGenes_manyAsTopTADs <- names(entrez_pval_sort[1:length(topTADs_genes)])
    
    stopifnot(length(topGenes_manyAsTopTADs) == length(topTADs_genes) )
    stopifnot(topTADs_genes %in% gene2tad_DT$entrezID)
    stopifnot(topTADs_genes %in% pipelineGenes)
    stopifnot(topGenes_manyAsTopPercent %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsTopPercent %in% pipelineGenes)
    stopifnot(topGenes_manyAsTopTADs %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsTopTADs %in% pipelineGenes)
    
    ### ENTROPY TEST
    # p_vals1 <- c(0.001,0.001,0.001,0.00005,0.001,0.002,0.002)
    # p_vals2 <- rep(0.00005, length(p_vals1))
    # 
    # p_vec1 <- -log10(p_vals1)
    # p_vec2 <- -log10(p_vals2)
    # 
    # proba_vec1 <- p_vec1/sum(p_vec1)
    # proba_vec2 <- p_vec2/sum(p_vec2)
    # 
    # proba_vec1 * -log2(proba_vec1) # -> use directly  proba_vec1 ???
    

        
    # run GO analysis for these 3 sets of entrez genes: topTADs_genes, topGenes_manyAsTopPercent, topGenes_manyAsTopTADs
    #***** 1) topTADs_genes
    topTADs_genes_enrich <- enricher(gene = topTADs_genes, 
                                               TERM2GENE=c5_msigdb,
                                               pvalueCutoff = enricher_pvalueCutoff, 
                                               pAdjustMethod = enricher_pAdjustMethod, 
                                               minGSSize = enricher_minGSSize, 
                                               maxGSSize = enricher_maxGSSize, 
                                               qvalueCutoff =enricher_qvalueCutoff)
    
    topTADs_genes_resultDT <- topTADs_genes_enrich@result
    topTADs_genes_resultDT <- topTADs_genes_resultDT[order(topTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    
    topTADs_genes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topTADs_genes_resultDT)))
    if(topTADs_genes_nTop > 0) {
      topTADs_genes_topGo  <- as.character(topTADs_genes_resultDT$ID[1:topTADs_genes_nTop])
    } else {
      topTADs_genes_topGo  <- character(0)  
    }
    stopifnot(length(topTADs_genes_topGo) == topTADs_genes_nTop)
    
    #***** 2) topGenes_manyAsTopPercent
    topGenes_manyAsTopPercent_enrich <- enricher(gene = topGenes_manyAsTopPercent, 
                                     TERM2GENE=c5_msigdb,
                                     pvalueCutoff = enricher_pvalueCutoff, 
                                     pAdjustMethod = enricher_pAdjustMethod, 
                                     minGSSize = enricher_minGSSize, 
                                     maxGSSize = enricher_maxGSSize, 
                                     qvalueCutoff =enricher_qvalueCutoff)
    
    topGenes_manyAsTopPercent_resultDT <- topGenes_manyAsTopPercent_enrich@result
    topGenes_manyAsTopPercent_resultDT <- topGenes_manyAsTopPercent_resultDT[order(topGenes_manyAsTopPercent_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    
    topGenes_manyAsTopPercent_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsTopPercent_resultDT)))
    if(topGenes_manyAsTopPercent_nTop > 0) {
      topGenes_manyAsTopPercent_topGo  <- as.character(topGenes_manyAsTopPercent_resultDT$ID[1:topGenes_manyAsTopPercent_nTop])
    } else {
      topGenes_manyAsTopPercent_topGo  <- character(0)  
    }
    stopifnot(length(topGenes_manyAsTopPercent_topGo) == topGenes_manyAsTopPercent_nTop)
    
    
    #***** 3) topGenes_manyAsTopTADs
    topGenes_manyAsTopTADs_enrich <- enricher(gene = topGenes_manyAsTopTADs, 
                                     TERM2GENE=c5_msigdb,
                                     pvalueCutoff = enricher_pvalueCutoff, 
                                     pAdjustMethod = enricher_pAdjustMethod, 
                                     minGSSize = enricher_minGSSize, 
                                     maxGSSize = enricher_maxGSSize, 
                                     qvalueCutoff =enricher_qvalueCutoff)

    topGenes_manyAsTopTADs_resultDT <- topGenes_manyAsTopTADs_enrich@result
    topGenes_manyAsTopTADs_resultDT <- topGenes_manyAsTopTADs_resultDT[order(topGenes_manyAsTopTADs_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    
    topGenes_manyAsTopTADs_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsTopTADs_resultDT)))
    if(topGenes_manyAsTopTADs_nTop > 0) {
      topGenes_manyAsTopTADs_topGo  <- as.character(topGenes_manyAsTopTADs_resultDT$ID[1:topGenes_manyAsTopTADs_nTop])
    } else {
      topGenes_manyAsTopTADs_topGo  <- character(0)  
    }
    stopifnot(length(topGenes_manyAsTopTADs_topGo) == topGenes_manyAsTopTADs_nTop)
    
    
    outFile <- file.path(outFold, paste0(curr_ds, "_enricher_all_results.Rdata"))
    enricher_all_results <- list(topTADs_genes = topTADs_genes_resultDT,
                                 topGenes_manyAsTopPercent = topGenes_manyAsTopPercent_resultDT,
                                 topGenes_manyAsTopTADs = topGenes_manyAsTopTADs_resultDT
                             )
    save(enricher_all_results, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    intersect_topTADs_topGenesManyAsTopPercent <- intersect(topGenes_manyAsTopPercent_topGo, topTADs_genes_topGo)
    intersect_topTADs_topGenesManyAsTopTADs <- intersect(topGenes_manyAsTopTADs_topGo, topTADs_genes_topGo)
    
    
    pvals0.01_topTADs_genes <- -log10(topTADs_genes_resultDT$p.adjust[topTADs_genes_resultDT$p.adjust <= 0.01])
    pvalsEntropy0.01_topTADs_genes <- mean(pvals0.01_topTADs_genes/sum(pvals0.01_topTADs_genes), na.rm=T)
    
    pvals0.05_topTADs_genes <- -log10(topTADs_genes_resultDT$p.adjust[topTADs_genes_resultDT$p.adjust <= 0.05])
    pvalsEntropy0.05_topTADs_genes <- mean(pvals0.05_topTADs_genes/sum(pvals0.05_topTADs_genes), na.rm=T)
    
    pvals0.01_topGenes_manyAsTopPercent <- -log10(topGenes_manyAsTopPercent_resultDT$p.adjust[topGenes_manyAsTopPercent_resultDT$p.adjust <= 0.01])
    pvalsEntropy0.01_topGenes_manyAsTopPercent <- mean(pvals0.01_topGenes_manyAsTopPercent/sum(pvals0.01_topGenes_manyAsTopPercent), na.rm=T)
    
    pvals0.05_topGenes_manyAsTopPercent <- -log10(topGenes_manyAsTopPercent_resultDT$p.adjust[topGenes_manyAsTopPercent_resultDT$p.adjust <= 0.05])
    pvalsEntropy0.05_topGenes_manyAsTopPercent <- mean(pvals0.05_topGenes_manyAsTopPercent/sum(pvals0.05_topGenes_manyAsTopPercent), na.rm=T)
    
    
    pvals0.01_topGenes_manyAsTopTADs <- -log10(topGenes_manyAsTopTADs_resultDT$p.adjust[topGenes_manyAsTopTADs_resultDT$p.adjust <= 0.01])
    pvalsEntropy0.01_topGenes_manyAsTopTADs <- mean(pvals0.01_topGenes_manyAsTopTADs/sum(pvals0.01_topGenes_manyAsTopTADs), na.rm=T)
    
    
    pvals0.05_topGenes_manyAsTopTADs <- -log10(topGenes_manyAsTopTADs_resultDT$p.adjust[topGenes_manyAsTopTADs_resultDT$p.adjust <= 0.05])
    pvalsEntropy0.05_topGenes_manyAsTopTADs <- mean(pvals0.05_topGenes_manyAsTopTADs/sum(pvals0.05_topGenes_manyAsTopTADs), na.rm=T)
    
    data.frame(
      dataset = curr_ds, 
      topTADs_genes_topGo = paste0(topTADs_genes_topGo, collapse=","),
      topGenes_manyAsTopPercent_topGo = paste0(topGenes_manyAsTopPercent_topGo, collapse=","),
      topGenes_manyAsTopTADs_topGo = paste0(topGenes_manyAsTopTADs_topGo, collapse=","),
      
      nTopPercentTADs = nTopTADs,
      nTopPercentTADs_genes = length(topTADs_genes),
      nTopPercentGenes = nTopGenes,
      
      nbrGO_topTADs_genes = topTADs_genes_nTop,
      nbrGO_topGenes_manyAsTopPercent =  topGenes_manyAsTopPercent_nTop,
      nbrGO_topGenes_manyAsTopTADs = topGenes_manyAsTopTADs_nTop,
      
      intersectNbrGO_topTADs_topGenesManyAsTopPercent = length(intersect_topTADs_topGenesManyAsTopPercent),
      intersectNbrGO_topTADs_topGenesManyAsTopTADs = length(intersect_topTADs_topGenesManyAsTopTADs),
      
      nResults_topTADs_genes = nrow(topTADs_genes_resultDT),
      nResults0.01_topTADs_genes = sum(topTADs_genes_resultDT$p.adjust <= 0.01),
      nResults0.05_topTADs_genes = sum(topTADs_genes_resultDT$p.adjust <= 0.05),
      
      nResults_topGenes_manyAsTopPercent = nrow(topGenes_manyAsTopPercent_resultDT),
      nResults0.01_topGenes_manyAsTopPercent = sum(topGenes_manyAsTopPercent_resultDT$p.adjust <= 0.01),
      nResults0.05_topGenes_manyAsTopPercent = sum(topGenes_manyAsTopPercent_resultDT$p.adjust <= 0.05),

      nResults_topGenes_manyAsTopTADs = nrow(topGenes_manyAsTopTADs_resultDT),
      nResults0.01_topGenes_manyAsTopTADs = sum(topGenes_manyAsTopTADs_resultDT$p.adjust <= 0.01),
      nResults0.05_topGenes_manyAsTopTADs = sum(topGenes_manyAsTopTADs_resultDT$p.adjust <= 0.05),
      
      pvalsEntropy0.01_topTADs_genes = pvalsEntropy0.01_topTADs_genes ,
      pvalsEntropy0.05_topTADs_genes = pvalsEntropy0.05_topTADs_genes,
      pvalsEntropy0.01_topGenes_manyAsTopPercent =  pvalsEntropy0.01_topGenes_manyAsTopPercent,
      pvalsEntropy0.05_topGenes_manyAsTopPercent = pvalsEntropy0.05_topGenes_manyAsTopPercent,
      pvalsEntropy0.01_topGenes_manyAsTopTADs = pvalsEntropy0.01_topGenes_manyAsTopTADs,
      pvalsEntropy0.05_topGenes_manyAsTopTADs =  pvalsEntropy0.05_topGenes_manyAsTopTADs,
      
      stringsAsFactors = FALSE
    )
    
    
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  orderDataset <- sort(setNames((all_ds_DT$intersectNbrGO_topTADs_topGenesManyAsTopPercent + all_ds_DT$intersectNbrGO_topTADs_topGenesManyAsTopTADs), all_ds_DT$dataset), decreasing=T)
  # orderDataset <- sort(setNames(pmax(all_ds_DT$intersect_topTADs_topGenesManyAsTopPercent, all_ds_DT$intersect_topTADs_topGenesManyAsTopTADs), all_ds_DT$dataset), decreasing=T)
  all_ds_DT <- all_ds_DT[match(names(orderDataset), all_ds_DT$dataset),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}


all_ds_DT$ratioTADgenes0.01 <- all_ds_DT$nResults0.01_topTADs_genes/all_ds_DT$nResults_topTADs_genes
all_ds_DT$ratioTADgenes0.05 <- all_ds_DT$nResults0.05_topTADs_genes/all_ds_DT$nResults_topTADs_genes
all_ds_DT$ratioPercentGenes0.01 <- all_ds_DT$nResults0.01_topGenes_manyAsTopPercent/all_ds_DT$nResults_topGenes_manyAsTopPercent
all_ds_DT$ratioPercentGenes0.05 <- all_ds_DT$nResults0.05_topGenes_manyAsTopPercent/all_ds_DT$nResults_topGenes_manyAsTopPercent
all_ds_DT$ratioManyTopTADsGenes0.01 <- all_ds_DT$nResults0.01_topGenes_manyAsTopTADs/all_ds_DT$nResults_topGenes_manyAsTopTADs
all_ds_DT$ratioManyTopTADsGenes0.05 <- all_ds_DT$nResults0.05_topGenes_manyAsTopTADs/all_ds_DT$nResults_topGenes_manyAsTopTADs

plot_multiDens(list(
  ratioTADgenes0.01 = all_ds_DT$ratioTADgenes0.01,
  ratioPercentGenes0.01 = all_ds_DT$ratioPercentGenes0.01,
  ratioManyTopTADsGenes0.01 = all_ds_DT$ratioManyTopTADsGenes0.01
),
plotTit = "Ratio signif. GO (p-val thresh = 0.01)"
)

plot_multiDens(list(
  ratioTADgenes0.05 = all_ds_DT$ratioTADgenes0.05,
  ratioPercentGenes0.05 = all_ds_DT$ratioPercentGenes0.05,
  ratioManyTopTADsGenes0.05 = all_ds_DT$ratioManyTopTADsGenes0.05
),
plotTit = "Ratio signif. GO (p-val thresh = 0.05)"
)

plot_multiDens(list(
  nResults_topTADs_genes = all_ds_DT$nResults_topTADs_genes,
  nResults_topGenes_manyAsTopPercent = all_ds_DT$nResults_topGenes_manyAsTopPercent,
  nResults_topGenes_manyAsTopTADs = all_ds_DT$nResults_topGenes_manyAsTopTADs
),
plotTit = "Tot. # enricher results",
  legPos = "topleft")


dsFold <- file.path(setDir, 
      "mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                    "OUTPUT_FOLDER")

all_ds_aucFCC <- foreach(curr_ds = all_ds_DT$dataset, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}
all_ds_DT$aucFCC <- all_ds_aucFCC

all_ds_aucCoexprDist <- foreach(curr_ds = all_ds_DT$dataset, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDist
}
all_ds_DT$aucCoexprDist <- all_ds_aucCoexprDist

dataset_order_FCC <- as.character(all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=T)])

my_plot_function <- function(var1, var2, DT, mylabels = NULL, withLab=F) {
  myx <- DT[,var1]
  myy <- DT[,var2]
  plot(x = myx,
       y = myy,
       xlab = ifelse(grepl("auc", var1), var1, paste0("signif. GO - ", var1)),
       ylab = paste0(var2),
       main = paste0(var2, " vs. ", var1),
       pch=16, cex=0.7
  )
  add_curv_fit(x=myx,
               y=myy, withR2 = F, lty=2)
  addCorr(x=myx,
          y=myy, bty="n")
  if(withLab){
    stopifnot(!is.null(mylabels))
    text(x = myx, y=myy, labels = mylabels, cex=0.6)
  }
}


all_curr_var <- c("ratioTADgenes0.05", "ratioTADgenes0.01", 
                  "ratioManyTopTADsGenes0.05", "ratioManyTopTADsGenes0.01", 
                  "ratioPercentGenes0.05", "ratioPercentGenes0.01",
                  "nResults_topTADs_genes",
                  "nResults_topGenes_manyAsTopTADs",
                  "nResults_topGenes_manyAsTopPercent",
                  "pvalsEntropy0.01_topTADs_genes",
                  "pvalsEntropy0.05_topTADs_genes",
                  "pvalsEntropy0.01_topGenes_manyAsTopPercent",
                  "pvalsEntropy0.05_topGenes_manyAsTopPercent",
                  "pvalsEntropy0.01_topGenes_manyAsTopTADs",
                  "pvalsEntropy0.05_topGenes_manyAsTopTADs"
                  )
all_ref_var <- c("aucFCC", "aucCoexprDist")

for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}




ref_var <- "ratioTADgenes0.05"
curr_var <- c("ratioManyTopTADsGenes0.05", "ratioPercentGenes0.05")
for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}



ref_var <- "ratioTADgenes0.01"
curr_var <- c("ratioManyTopTADsGenes0.01", "ratioPercentGenes0.01")
for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}


ref_var <- "nResults_topTADs_genes"
curr_var <- c("nResults_topGenes_manyAsTopTADs", "nResults_topGenes_manyAsTopPercent")
for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifGO_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}



########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES

all_pval <- c("0.01", "0.05")

all_types <- c("PercentGenes", "ManyTopTADsGenes")

pval=0.01
curr_type="ManyTopTADs"

for(pval in all_pval){
  for(curr_type in all_types){
    # if(curr_type == )
    plotDT <- all_ds_DT[,c("dataset", paste0("ratioTADgenes", pval), paste0("ratio",  curr_type, pval))]
    plotDT_m <- melt(plotDT, id="dataset")
    
    datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
    plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
    stopifnot(!is.na(plotDT_m$dataset))
    
    plotTit <- paste0("Ratio signif GOs ", pval)
    mySub <- paste0("top Genes = ", curr_type)
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      
      scale_fill_manual(values =setNames( c("dodgerblue4", "darkorange2"), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ),
                        labels =setNames( c("TADgenes", curr_type), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ) ) +
      labs(fill = "")+
      theme( # Increase size of axis lines
        # top, right, bottom and left
        plot.margin = unit(c(1, 1, 2, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=8),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank()
      )
    if(SSHFS) p_var
    outFile <- file.path(outFold, paste0("signifGO_ratioTADgenes", pval, "_ratio", curr_type, pval, ".", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}
########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES - order FCC

all_pval <- c("0.01", "0.05")

all_types <- c("PercentGenes", "ManyTopTADsGenes")

pval=0.01
curr_type="ManyTopTADs"

for(pval in all_pval){
  for(curr_type in all_types){
    # if(curr_type == )
    plotDT <- all_ds_DT[,c("dataset", paste0("ratioTADgenes", pval), paste0("ratio",  curr_type, pval))]
    plotDT_m <- melt(plotDT, id="dataset")
    
    # datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
    datasetOrder <- dataset_order_FCC
    plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
    stopifnot(!is.na(plotDT_m$dataset))
    
    plotTit <- paste0("Ratio signif GOs ", pval)
    mySub <- paste0("top Genes = ", curr_type)
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      
      scale_fill_manual(values =setNames( c("dodgerblue4", "darkorange2"), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ),
                        labels =setNames( c("TADgenes", curr_type), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ) ) +
      labs(fill = "")+
      theme( # Increase size of axis lines
        # top, right, bottom and left
        plot.margin = unit(c(1, 1, 2, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=8),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank()
      )
    if(SSHFS) p_var
    outFile <- file.path(outFold, paste0("signifGO_ratioTADgenes", pval, "_ratio", curr_type, pval, "_orderFCC.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}


########################################## BARPLOT ENTROPY BOTH KIND OF # OF GENES - order FCC

# pvalsEntropy0.01_topGenes_manyAsTopPercent
all_pval <- c("0.01", "0.05")

all_types <- c("topGenes_manyAsTopPercent", "topGenes_manyAsTopTADs")

pval=0.01
curr_type="ManyTopTADs"

for(pval in all_pval){
  for(curr_type in all_types){
    # if(curr_type == )
    
    var1 <- paste0("pvalsEntropy", pval, "_topTADs_genes")
    var2 <- paste0("pvalsEntropy", pval,"_", curr_type)
    
    plotDT <- all_ds_DT[,c("dataset", var1, var2)]
    plotDT_m <- melt(plotDT, id="dataset")
    
    # datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
    datasetOrder <- dataset_order_FCC
    plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
    stopifnot(!is.na(plotDT_m$dataset))
    
    plotTit <- paste0("Mean pvals entropy ", pval)
    
    mySub <- paste0("top Genes = ", curr_type)
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values =  setNames(c("dodgerblue4", "darkorange2"), c(var1, var2)),
                        labels =  setNames(c("topTADs_genes", curr_type), c(var1, var2)))+
      # scale_fill_manual(values =setNames( c("dodgerblue4", "darkorange2"), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ),
      #                   labels =setNames( c("TADgenes", curr_type), c(paste0("ratioTADgenes", pval), paste0("ratio", curr_type, pval)) ) ) +
      labs(fill = "")+
      theme( # Increase size of axis lines
        # top, right, bottom and left
        plot.margin = unit(c(1, 1, 2, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=8),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank()
      )
    if(SSHFS) p_var
    outFile <- file.path(outFold, paste0("signifGO_pvalsEntropyTADgenes", pval, "_pvalsEntropy", curr_type, pval, "_orderFCC.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}




######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
