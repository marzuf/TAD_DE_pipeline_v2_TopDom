startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO_pvalSelect.R\n"))

# Rscript cmp_DE_TADs_genes_GO_pvalSelect.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

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

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05

outFold <- "CMP_TADs_GENES_GO_pvalSelect"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

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

#all_ds <- all_ds[1:3]

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
    

    stopifnot(selectTADs %in% gene2tad_DT$region)
    
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
    
    
    
    # run GO analysis for these 3 sets of entrez genes: selectTADs_genes, selectGenes
    #***** 1) selectTADs_genes
    if(nSelectTADs > 0) {
      selectTADs_genes_enrich <- enricher(gene = selectTADs_genes, 
                                          TERM2GENE=c5_msigdb,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
      
      selectTADs_genes_resultDT <- selectTADs_genes_enrich@result
      selectTADs_genes_resultDT <- selectTADs_genes_resultDT[order(selectTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
      
      selectTADs_genes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(selectTADs_genes_resultDT)))
      if(selectTADs_genes_nTop > 0) {
        selectTADs_genes_topGo  <- as.character(selectTADs_genes_resultDT$ID[1:selectTADs_genes_nTop])
      } else {
        selectTADs_genes_topGo  <- character(0)  
      }
      stopifnot(length(selectTADs_genes_topGo) == selectTADs_genes_nTop)
    } else {
      selectTADs_genes_resultDT <- data.frame(p.adjust = NA)
      selectTADs_genes_topGo <- NA
      selectTADs_genes_nTop <- NA
      
    }

    
    #***** 2) selectGenes
    if(nSelectGenes > 0) {
      selectGenes_enrich <- enricher(gene = selectGenes, 
                                     TERM2GENE=c5_msigdb,
                                     pvalueCutoff = enricher_pvalueCutoff, 
                                     pAdjustMethod = enricher_pAdjustMethod, 
                                     minGSSize = enricher_minGSSize, 
                                     maxGSSize = enricher_maxGSSize, 
                                     qvalueCutoff =enricher_qvalueCutoff)
      
      selectGenes_resultDT <- selectGenes_enrich@result
      selectGenes_resultDT <- selectGenes_resultDT[order(selectGenes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
      
      selectGenes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(selectGenes_resultDT)))
      if(selectGenes_nTop > 0) {
        selectGenes_topGo  <- as.character(selectGenes_resultDT$ID[1:selectGenes_nTop])
      } else {
        selectGenes_topGo  <- character(0)  
      }
      stopifnot(length(selectGenes_topGo) == selectGenes_nTop)
    } else {
      
      selectGenes_resultDT <- data.frame(p.adjust = NA)
      selectGenes_topGo <- NA
      selectGenes_nTop <- NA
    }

    
    
    outFile <- file.path(outFold, paste0(curr_ds, "_enricher_all_results.Rdata"))
    enricher_all_results <- list(selectTADs_genes = selectTADs_genes_resultDT,
                                 selectGenes = selectGenes_resultDT
                             )
    save(enricher_all_results, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    intersect_selectTADs_genes_selectGenes <- intersect(selectGenes_topGo, selectTADs_genes_topGo)
    
    data.frame(
      dataset = curr_ds, 
      selectTADs_genes_topGo = paste0(selectTADs_genes_topGo, collapse=","),
      selectGenes_topGo = paste0(selectGenes_topGo, collapse=","),
      
      nSelectTADs = nSelectTADs,
      nSelectTADs_genes = length(selectTADs_genes),
      nSelectGenes = nSelectGenes,
      
      nbrGO_selectTADs_genes = selectTADs_genes_nTop,
      nbrGO_selectGenes =  selectGenes_nTop,
      
      intersectNbrGO_selectTADs_genes_selectGenes = length(intersect_selectTADs_genes_selectGenes),
      
      stringsAsFactors = FALSE
    )
    
    
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  orderDataset <- sort(setNames((all_ds_DT$intersectNbrGO_selectTADs_genes_selectGenes), all_ds_DT$dataset), decreasing=T)

  all_ds_DT <- all_ds_DT[match(names(orderDataset), all_ds_DT$dataset),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}


all_ds_DT$nGO_selectTADs_genes_selectGenes <- (all_ds_DT$nbrGO_selectTADs_genes + all_ds_DT$nbrGO_selectGenes)
all_ds_DT$ratio_selectGenes <- all_ds_DT$intersectNbrGO_selectTADs_genes_selectGenes/all_ds_DT$nGO_selectTADs_genes_selectGenes

stopifnot(na.omit(all_ds_DT$ratio_selectGenes) >= 0 & na.omit(all_ds_DT$ratio_selectGenes) <= 1)


########################################## COMPARE % TOP GENES AND SAME # GENES AS TOP TADS
myx <- all_ds_DT$nSelectTADs_genes
myy <- all_ds_DT$nSelectGenes
outFile <- file.path(outFold, paste0("cmp_totNbr_diffNbrGenes.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "# genes compared",
     xlab="genes from DE TADs",
     ylab="DE genes",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES

plotDT <- all_ds_DT[,c("dataset", "ratio_selectGenes")]
plotDT_m <- melt(plotDT, id="dataset")

datasetOrder <- plotDT$dataset[order(plotDT[,2],  decreasing=T)]
plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
stopifnot(!is.na(plotDT_m$dataset))

plotTit <- "Ratio intersect GOs"

p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
  ggtitle(plotTit, subtitle = mySub)+
  geom_bar(stat="identity", position = "dodge")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(""),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values = c(ratio_selectGenes = "dodgerblue4"),
                    labels = c(ratio_selectGenes = "DE genes"))+
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

outFile <- file.path(outFold, paste0("GO_at_intersect_", "bothRatios.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

########################################## BARPLOT # VS. # INTERSECT FOR EACH KIND OF # OF GENES

toplot <- ""

for(curr_plot in toplot) {
    plotDT <- all_ds_DT[,c("dataset", "nGO_selectTADs_genes_selectGenes", "intersectNbrGO_selectTADs_genes_selectGenes")]
    plotTit <- paste0("# GO and intersect - select DE genes")
    
  plotDT_m <- melt(plotDT, id="dataset")
  
  datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
  plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
  stopifnot(!is.na(plotDT_m$dataset))
  
  p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
    ggtitle(plotTit, subtitle = mySub)+
    geom_bar(stat="identity", position = "dodge")+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 10))+
    guides(fill=FALSE)+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
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
  
  outFile <- file.path(outFold, paste0("GO_at_intersect_", curr_plot, "TopGenes.", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
