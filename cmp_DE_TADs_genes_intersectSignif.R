startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_intersectSignif.R\n"))

# Rscript cmp_DE_TADs_genes_intersectSignif.R


options(scipen=100)

library(foreach)
library(doMC)
library(reshape2)
library(ggplot2)

source("analysis_utils.R")

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

buildTable <- FALSE

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

outFold <- "CMP_TADs_GENES_INTERSECT"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_intersect_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

# all_ds <- all_ds[1:10]

if(buildTable) {
  
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
    
    intersect_manyAsPercent <- intersect(topTADs_genes, topGenes_manyAsTopPercent)
    intersect_manyAsTopTADs <- intersect(topTADs_genes, topGenes_manyAsTopTADs)
    
    data.frame(
      dataset = curr_ds, 
      nTopPercentTADs = nTopTADs,
      nTopPercentTADs_genes = length(topTADs_genes),
      nTopPercentGenes = nTopGenes,
      intersectNbrGenes_nbrPercent = length(intersect_manyAsPercent),
      intersectNbrGenes_sameNbr = length(intersect_manyAsTopTADs),
      stringsAsFactors = FALSE
    )
  } # end iterate over all datasets
  rownames(all_ds_DT) <- NULL
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

all_ds_DT$ratio_nbrPercent <- all_ds_DT$intersectNbrGenes_nbrPercent/all_ds_DT$nTopPercentGenes
all_ds_DT$ratio_sameNbr <- all_ds_DT$intersectNbrGenes_sameNbr/all_ds_DT$nTopPercentTADs_genes
stopifnot(all_ds_DT$ratio_nbrPercent >= 0 & all_ds_DT$ratio_nbrPercent <= 1)
stopifnot(all_ds_DT$ratio_sameNbr >= 0 & all_ds_DT$ratio_sameNbr <= 1)

########################################## COMPARE % TOP GENES AND SAME # GENES AS TOP TADS

myx <- all_ds_DT$nTopPercentTADs_genes
myy <- all_ds_DT$nTopPercentGenes
outFile <- file.path(outFold, paste0("cmp_totNbr_diffNbrGenes.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "# genes compared",
     xlab="as many genes as from topTADs",
     ylab=paste0("top ", round(topPercent * 100, 2), "% DE genes"),
     pch=16,cex=0.7)
add_curv_fit(x=myx,
    y=myy, withR2 = F, lty=2)
addCorr(x=myx,
          y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

    
myx <- all_ds_DT$intersectNbrGenes_sameNbr
myy <- all_ds_DT$intersectNbrGenes_nbrPercent
outFile <- file.path(outFold, paste0("cmp_intersect_diffNbrGenes.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "# genes intersect",
     xlab="intersect as many genes as from topTADs",
     ylab=paste0("intersect top ", round(topPercent * 100, 2), "% DE genes"),
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES

plotDT <- all_ds_DT[,c("dataset", "ratio_nbrPercent", "ratio_sameNbr")]
plotDT_m <- melt(plotDT, id="dataset")

datasetOrder <- plotDT$dataset[order(plotDT[,2], plotDT[,3], decreasing=T)]
plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
stopifnot(!is.na(plotDT_m$dataset))

plotTit <- "Ratio intersect genes"

p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
  ggtitle(plotTit)+
  geom_bar(stat="identity", position = "dodge")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(""),
                     breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = c(ratio_nbrPercent = "dodgerblue4", ratio_sameNbr = "darkorange2"),
                    labels = c(ratio_nbrPercent = "% top genes", ratio_sameNbr = "same # genes as topTADs"))+
  labs(fill = "")+
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

outFile <- file.path(outFold, paste0("genes_at_intersect_bothRatios.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

########################################## BARPLOT # VS. # INTERSECT FOR EACH KIND OF # OF GENES

toplot <- c("percent", "sameNbr") 

for(curr_plot in toplot) {
  if(curr_plot=="percent") {
    plotDT <- all_ds_DT[,c("dataset", "nTopPercentGenes", "intersectNbrGenes_nbrPercent")]
    plotTit <- paste0("# genes and intersect - top ", round(topPercent*100, 2), "% DE genes")
  } else if(curr_plot=="sameNbr") {
    plotDT <- all_ds_DT[,c("dataset", "nTopPercentTADs_genes", "intersectNbrGenes_sameNbr")]
    plotTit <- paste0("# genes and intersect - as many genes as in top TADs")
  }
  plotDT_m <- melt(plotDT, id="dataset")

  datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
  plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
  stopifnot(!is.na(plotDT_m$dataset))

    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit)+
    geom_bar(stat="identity", position = "dodge")+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 5))+
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
    outFile <- file.path(outFold, paste0("genes_at_intersect_", curr_plot, "TopGenes.", plotType))
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
