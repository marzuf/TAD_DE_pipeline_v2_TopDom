startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_sum_ranks.R\n"))
# Rscript cmp_DE_TADs_genes_sum_ranks.R

options(scipen=100)

library(foreach)
library(doMC)
library(reshape2)
library(ggplot2)

buildTable <- F

plotType <- "png"
myHeight <- ifelse(plotType == "png", 500, 7)
myWidth <- myHeight


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source("analysis_utils.R")

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

outFold <- "CMP_SUM_RANKS"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_sum_ranks_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

toPrintTADs <- c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")
toPrintDs <- c("TCGAcrc_msi_mss")

# all_ds <- all_ds[1:10]

#**************************** hard-coded settings
### PLOT TOP 5:
dataset_to_plot <- c("TCGAcrc_msi_mss", "GSE102073_stic_nostic", "GSE74927_neg_pos", "TCGAacc_acc_mutCTNNB1", "TCGAskcm_skcm_mutCTNNB1")

#**************************************************



if(buildTable){
  all_ds_all_TADs_DT <- foreach(curr_ds = all_ds, .combine='rbind') %do% {
    
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
    
    # for each TAD
    # -> TAD rank
    # -> genes that are in TAD -> sum gene rank
    
    all_tads <- names(tad_pval_rank)
    
    all_TADs_DT <- foreach(curr_tad = all_tads, .combine='rbind') %dopar% {
      
      TAD_rank <- tad_pval_rank[curr_tad]
      
      TAD_genes <- as.character(gene2tad_DT$entrezID[gene2tad_DT$region == curr_tad])
      TAD_genes <- TAD_genes[TAD_genes %in% names(entrez_pval_rank)]
      
      stopifnot(length(TAD_genes) > 0)
      
      geneRanks <- entrez_pval_rank[as.character(TAD_genes)]
      geneRanks <- sort(geneRanks)
      
      stopifnot(!is.na(geneRanks))
      stopifnot(!is.na(names(geneRanks)))
      
      data.frame(
        TAD = curr_tad,
        TAD_rank = TAD_rank,
        TAD_genes = paste0(names(geneRanks), collapse=","),
        TAD_genes_rank = paste0(geneRanks, collapse=","),
        nGenes = length(geneRanks),
        sumGeneRank = sum(geneRanks),
        stringsAsFactors = FALSE
      )
      
    } # end iterate over all TADs
    all_TADs_DT$dataset <- curr_ds
    all_TADs_DT <- all_TADs_DT[, c("dataset", colnames(all_TADs_DT)[colnames(all_TADs_DT) != "dataset"])]
    all_TADs_DT$genes_meanRank <- all_TADs_DT$sumGeneRank/all_TADs_DT$nGenes
    all_TADs_DT
  } # end iterate over all datasets
  
  rownames(all_ds_all_TADs_DT) <- NULL
  
  outFile <- file.path(outFold, "all_ds_all_TADs_DT.Rdata")
  save(all_ds_all_TADs_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_all_TADs_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_all_TADs_DT <- eval(parse(text = load(outFile)))
}

#### ALL DATASETS  ################################################################
myx <- all_ds_all_TADs_DT$TAD_rank
myy <- all_ds_all_TADs_DT$genes_meanRank
outFile <- file.path(outFold, paste0("cmp_TADrank_genesMeanRank.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "Rank TAD vs. rank TAD genes",
     xlab="TAD rank",
     ylab="mean rank TAD genes",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
mtext(text = paste0("all datasets (n=", length(unique(all_ds_all_TADs_DT$dataset)), "), all TADs (n=",nrow(all_ds_all_TADs_DT), ")") , 
      side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- all_ds_all_TADs_DT$TAD_rank
myy <- all_ds_all_TADs_DT$sumGeneRank
outFile <- file.path(outFold, paste0("cmp_TADrank_genesRankSum.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "Rank TAD vs. sum TAD genes rank",
     xlab="TAD rank",
     ylab="sum TAD genes rank",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
mtext(text = paste0("all datasets (n=", length(unique(all_ds_all_TADs_DT$dataset)), "), all TADs (n=",nrow(all_ds_all_TADs_DT), ")") , 
      side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

### dataset_to_plot ################################################################
toplotDT <- all_ds_all_TADs_DT[all_ds_all_TADs_DT$dataset %in% dataset_to_plot,]
myx <- toplotDT$TAD_rank
myy <- toplotDT$genes_meanRank
outFile <- file.path(outFold, paste0("cmp_TADrank_genesMeanRank_datasetToPlot.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "Rank TAD vs. rank TAD genes",
     xlab="TAD rank",
     ylab="mean rank TAD genes",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
mtext(text = paste0(paste0(dataset_to_plot, collapse = ", "), " - all TADs (n=",nrow(all_ds_all_TADs_DT), ")") , 
      side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- toplotDT$TAD_rank
myy <- toplotDT$sumGeneRank
outFile <- file.path(outFold, paste0("cmp_TADrank_genesRankSum_datasetToPlot.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "Rank TAD vs. sum TAD genes rank",
     xlab="TAD rank",
     ylab="sum TAD genes rank",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
mtext(text = paste0(paste0(dataset_to_plot, collapse = ", "), " - all TADs (n=",nrow(all_ds_all_TADs_DT), ")") , 
      side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


toPrint_DT <- all_ds_all_TADs_DT[all_ds_all_TADs_DT$dataset %in% toPrintDs & all_ds_all_TADs_DT$TAD %in% toPrintTADs,]
printAndLog("\n", logFile)
write.table(toPrint_DT, file = logFile, col.names=T, row.names=F, append=T, sep="\t", quote=F)

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
