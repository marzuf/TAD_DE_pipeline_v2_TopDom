startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_KEGGsignif_pvalSelect.R\n"))

# Rscript cmp_DE_TADs_genes_KEGGsignif_pvalSelect.R

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

# topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05

outFold <- "CMP_TADs_GENES_KEGG_SIGNIF_pvalSelect"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_kegg_signif_pvalSelect_logFile.txt")
system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)


# settings for enrichment analysis: clusterProfiler::enrichKEGG
# default values enrichKEGG function: (universe as background !)
# pvalueCutoff = 0.05 
# pAdjustMethod = "BH"
# minGSSize = 10
# maxGSSize = 500
# qvalueCutoff = 0.2 
enrichKEGG_pvalueCutoff <- 1
enrichKEGG_pAdjustMethod <- "BH"
enrichKEGG_minGSSize <- 1
enrichKEGG_maxGSSize <- 500
enrichKEGG_qvalueCutoff <- 1
# logical, use KEGG.db or latest online KEGG data
enrichKEGG_internal <- FALSE

enrichKEGG_results_sortKEGGby <- "p.adjust"
enrichKEGG_results_cmpNbrTop <- 50

mySub <- paste0("(# top KEGG: ",enrichKEGG_results_cmpNbrTop, ")") 

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG pvalueCutoff:\t", enrichKEGG_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG pAdjustMethod:\t", enrichKEGG_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG minGSSize:\t", enrichKEGG_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG maxGSSize:\t", enrichKEGG_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG qvalueCutoff:\t", enrichKEGG_qvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enrichKEGG use_internal_data:\t", enrichKEGG_internal, "\n")
printAndLog(txt, logFile)

#all_ds=all_ds[1:3]
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
    
    # -> change here pval selection
    selectTADs <- names(tad_pval[tad_pval <= pvalSelect])
    nSelectTADs <- length(selectTADs)
    
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
    selectGenes <- names(entrez_pval[entrez_pval <= pvalSelect])
    stopifnot(selectGenes %in% pipelineGenes)
    nSelectGenes <- length(selectGenes)

    if(nSelectTADs > 0){
      stopifnot(selectTADs %in% gene2tad_DT$region)
      selectTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% selectTADs]  
      selectTADs_genes <- selectTADs_genes[selectTADs_genes %in% pipelineGenes]
      stopifnot(selectTADs_genes %in% gene2tad_DT$entrezID)
      stopifnot(selectTADs_genes %in% pipelineGenes)
      
    }
  
  if(nSelectGenes > 0){
    stopifnot(selectGenes %in% gene2tad_DT$entrezID)
    stopifnot(selectGenes %in% pipelineGenes)
  }
  
    # run KEGG analysis for these 2 sets of entrez genes: selectTADs_genes, selectGenes,
    #***** 1) selectTADs_genes
    if(nSelectTADs > 0){
      
      selectTADs_genes_enrich <- enrichKEGG(gene = selectTADs_genes, 
                                         organism = "hsa", keyType = "kegg",
                                         pvalueCutoff = enrichKEGG_pvalueCutoff, 
                                         pAdjustMethod = enrichKEGG_pAdjustMethod, 
                                         minGSSize = enrichKEGG_minGSSize, 
                                         maxGSSize = enrichKEGG_maxGSSize, 
                                         qvalueCutoff =enrichKEGG_qvalueCutoff,
                                         use_internal_data = enrichKEGG_internal
                                          )
      
      selectTADs_genes_resultDT <- selectTADs_genes_enrich@result
      selectTADs_genes_resultDT <- selectTADs_genes_resultDT[order(selectTADs_genes_resultDT[,enrichKEGG_results_sortKEGGby], decreasing=FALSE), ]
      
      selectTADs_genes_nTop <- min(c(enrichKEGG_results_cmpNbrTop, nrow(selectTADs_genes_resultDT)))
      if(selectTADs_genes_nTop > 0) {
        selectTADs_genes_topKEGG  <- as.character(selectTADs_genes_resultDT$ID[1:selectTADs_genes_nTop])
      } else {
        selectTADs_genes_topKEGG  <- character(0)  
      }
      stopifnot(length(selectTADs_genes_topKEGG) == selectTADs_genes_nTop)
    } else {
      selectTADs_genes_resultDT <- data.frame(p.adjust = NA)
      selectTADs_genes_topKEGG <- NA
      selectTADs_genes_nTop <- NA
    }
    #***** 2) selectGenes
    if(nSelectGenes > 0){
      selectGenes_enrich <- enrichKEGG(gene = selectGenes, 
                                                     organism = "hsa", keyType = "kegg",
                                                     pvalueCutoff = enrichKEGG_pvalueCutoff, 
                                                     pAdjustMethod = enrichKEGG_pAdjustMethod, 
                                                     minGSSize = enrichKEGG_minGSSize, 
                                                     maxGSSize = enrichKEGG_maxGSSize, 
                                                     qvalueCutoff =enrichKEGG_qvalueCutoff,
                                                     use_internal_data = enrichKEGG_internal)
      
      selectGenes_resultDT <- selectGenes_enrich@result
      selectGenes_resultDT <- selectGenes_resultDT[order(selectGenes_resultDT[,enrichKEGG_results_sortKEGGby], decreasing=FALSE), ]
      
      selectGenes_nTop <- min(c(enrichKEGG_results_cmpNbrTop, nrow(selectGenes_resultDT)))
      if(selectGenes_nTop > 0) {
        selectGenes_topKEGG  <- as.character(selectGenes_resultDT$ID[1:selectGenes_nTop])
      } else {
        selectGenes_topKEGG  <- character(0)  
      }
      stopifnot(length(selectGenes_topKEGG) == selectGenes_nTop)
    } else {
      selectGenes_resultDT <- data.frame(p.adjust = NA)
      selectGenes_topKEGG <- NA
      selectGenes_nTop <- NA
    }
    
    outFile <- file.path(outFold, paste0(curr_ds, "_enrichKEGG_all_results.Rdata"))
    enrichKEGG_all_results <- list(selectTADs_genes = selectTADs_genes_resultDT,
                                 selectGenes = selectGenes_resultDT)
    
    save(enrichKEGG_all_results, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    if(nSelectGenes == 0 | nSelectTADs == 0){
      intersect_selectTADs_genes_selectGenes <- NA
    } else {
      intersect_selectTADs_genes_selectGenes <- intersect(selectGenes_topKEGG, selectTADs_genes_topKEGG)  
    }
    
    
    data.frame(
      dataset = curr_ds, 
      selectTADs_genes_topKEGG = paste0(selectTADs_genes_topKEGG, collapse=","),
      selectGenes_topKEGG = paste0(selectGenes_topKEGG, collapse=","),
      
      nSelectTADs = nSelectTADs,
      nSelectGenes = nSelectGenes,
      
      nbrKEGG_selectTADs_genes = selectTADs_genes_nTop,
      nbrKEGG_selectGenes =  selectGenes_nTop,
      
      intersectNbrKEGG_selectTADs_genes_selectGenes = length(intersect_selectTADs_genes_selectGenes),
      
      nResults_selectTADs_genes = nrow(selectTADs_genes_resultDT),
      nResults0.01_selectTADs_genes = sum(selectTADs_genes_resultDT$p.adjust <= 0.01),
      nResults0.05_selectTADs_genes = sum(selectTADs_genes_resultDT$p.adjust <= 0.05),
      
      nResults_selectGenes = nrow(selectGenes_resultDT),
      nResults0.01_selectGenes = sum(selectGenes_resultDT$p.adjust <= 0.01),
      nResults0.05_selectGenes = sum(selectGenes_resultDT$p.adjust <= 0.05),
      
      stringsAsFactors = FALSE
    )
    
    
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  orderDataset <- sort(setNames((all_ds_DT$intersectNbrKEGG_selectTADs_genes_selectGenes), all_ds_DT$dataset), decreasing=T)
  all_ds_DT <- all_ds_DT[match(names(orderDataset), all_ds_DT$dataset),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFold <- "CMP_TADs_GENES_KEGG_SIGNIF_pvalSelect"
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

all_ds_DT$ratioTADgenes0.01 <- all_ds_DT$nResults0.01_selectTADs_genes/all_ds_DT$nResults_selectTADs_genes
all_ds_DT$ratioTADgenes0.05 <- all_ds_DT$nResults0.05_selectTADs_genes/all_ds_DT$nResults_selectTADs_genes
all_ds_DT$ratioSelectGenes0.01 <- all_ds_DT$nResults0.01_selectGenes/all_ds_DT$nResults_selectGenes
all_ds_DT$ratioSelectGenes0.05 <- all_ds_DT$nResults0.05_selectGenes/all_ds_DT$nResults_selectGenes

stopifnot(na.omit(all_ds_DT$ratioTADgenes0.01) >= 0 & na.omit(all_ds_DT$ratioTADgenes0.01) <= 1)
stopifnot(na.omit(all_ds_DT$ratioTADgenes0.05) >= 0 & na.omit(all_ds_DT$ratioTADgenes0.05) <= 1)
stopifnot(na.omit(all_ds_DT$ratioSelectGenes0.01) >= 0 & na.omit(all_ds_DT$ratioSelectGenes0.01) <= 1)
stopifnot(na.omit(all_ds_DT$ratioSelectGenes0.05) >= 0 & na.omit(all_ds_DT$ratioSelectGenes0.05) <= 1)


plot_multiDens(list(
  ratioTADgenes0.01 = all_ds_DT$ratioTADgenes0.01,
  ratioSelectGenes0.01 = all_ds_DT$ratioSelectGenes0.01
),
plotTit = "Ratio signif. KEGG (p-val thresh = 0.01)"
)

plot_multiDens(list(
  ratioTADgenes0.05 = all_ds_DT$ratioTADgenes0.05,
  ratioSelectGenes0.05 = all_ds_DT$ratioSelectGenes0.05
),
plotTit = "Ratio signif. KEGG (p-val thresh = 0.05)"
)

plot_multiDens(list(
  nResults_selectTADs_genes = all_ds_DT$nResults_selectTADs_genes,
  nResults_selectGenes = all_ds_DT$nResults_selectGenes
  ),
plotTit = "Tot. # enrichKEGG results",
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
       xlab = ifelse(grepl("auc", var1), var1, paste0("signif. KEGG - ", var1)),
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
                  "ratioSelectGenes0.05", "ratioSelectGenes0.01",
                  "nResults_selectTADs_genes",
                  "nResults_selectGenes"
                  )
stopifnot(all_curr_var %in% colnames(all_ds_DT))

all_ref_var <- c("aucFCC", "aucCoexprDist")
stopifnot(all_ref_var %in% colnames(all_ds_DT))

for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifKEGG.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifKEGG_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}


ref_var <- "nResults_selectTADs_genes"
stopifnot(all_curr_var %in% colnames(all_ds_DT))

curr_var <- c("nResults_selectGenes")
stopifnot(curr_var %in% colnames(all_ds_DT))

for(ref_var in all_ref_var){
  for(curr_var in all_curr_var) {
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifKEGG.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, withLab=FALSE)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(ref_var, "_vs_", curr_var, "_signifKEGG_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(ref_var, curr_var, all_ds_DT, mylabels = as.character(all_ds_DT$dataset), withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}



########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES

all_pval <- c("0.01", "0.05")

all_types <- c("SelectGenes")

pval=0.01
curr_type="SelectGenes"

for(pval in all_pval){
  for(curr_type in all_types){
    # if(curr_type == )
    plotDT <- all_ds_DT[,c("dataset", paste0("ratioTADgenes", pval), paste0("ratio",  curr_type, pval))]
    plotDT_m <- melt(plotDT, id="dataset")
    
    datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
    plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
    stopifnot(!is.na(plotDT_m$dataset))
    
    plotTit <- paste0("Ratio signif KEGGs ", pval)
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
    outFile <- file.path(outFold, paste0("signifKEGG_ratioTADgenes", pval, "_ratio", curr_type, pval, ".", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}
########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES - order FCC

all_pval <- c("0.01", "0.05")

all_types <- c("SelectGenes")

pval=0.01
curr_type="ManyTopTADs"

for(pval in all_pval){
  for(curr_type in all_types){
    # if(curr_type == )
    plotDT <- all_ds_DT[,c("dataset", paste0("ratioTADgenes", pval), paste0("ratio",  curr_type, pval))]
    plotDT_m <- melt(plotDT, id="dataset")
    
    # datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,1], decreasing=T)]
    datasetOrder <- dataset_order_FCC
    plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
    stopifnot(!is.na(plotDT_m$dataset))
    
    plotTit <- paste0("Ratio signif KEGGs ", pval)
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
    outFile <- file.path(outFold, paste0("signifKEGG_ratioTADgenes", pval, "_ratio", curr_type, pval, "_orderFCC.", plotType))
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
