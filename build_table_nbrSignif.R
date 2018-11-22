startTime <- Sys.time()
cat(paste0("> build_table_nbrSignif.R\n"))

# Rscript build_table_nbrSignif.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight


buildTable <- TRUE

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

source("analysis_utils.R")


outFold <- "BUILD_TABLE_NBR_SIGNIF"
system(paste0("mkdir -p ", outFold))

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

all_pval_thresh <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
sub_pval_thresh <- c(0.01, 0.05)

if(buildTable) {
  all_ds_DT <- foreach(curr_ds = all_ds, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    
    ### RETRIEVE TAD PVALUES
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    tad_pval <- sort(tad_pval, decreasing=F)
    stopifnot(is.numeric(tad_pval))
    nSignifTADs <- sapply(all_pval_thresh, function(pval) sum(tad_pval <= pval))
    nTADs <- length(tad_pval)
    
    ### RETRIEVE GENES PVALUES
    step1_fold <- file.path(dsFold, curr_ds, "1_runGeneDE")
    stopifnot(file.exists(step1_fold))
    toptableFile <- file.path(step1_fold, "DE_topTable.Rdata")
    stopifnot(file.exists(toptableFile))
    topTable_DT <- eval(parse(text = load(toptableFile)))
    stopifnot(!any(duplicated(topTable_DT$genes)))
    stopifnot(topTable_DT$genes == rownames(topTable_DT))
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(names(pipelineGenes) %in% topTable_DT$genes)
    topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipelineGenes),]
    stopifnot(nrow(topTable_DT) > 0)
    stopifnot(is.numeric(topTable_DT$adj.P.Val))
    nSignifGenes <- sapply(all_pval_thresh, function(pval) sum(topTable_DT$adj.P.Val <= pval))
    nGenes <- nrow(topTable_DT)
    
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
    
    all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
    aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
    stopifnot(!is.na(aucCoexprDist))
    
    stopifnot(length(nSignifGenes) == length(nSignifTADs))
    stopifnot(length(nSignifGenes) == length(all_pval_thresh))
    
    c(
      nTADs,
      nSignifTADs,
      nGenes,
      nSignifGenes,
      aucFCC,
      aucCoexprDist
    )  
  }
  rownames(all_ds_DT) <- all_ds
  all_ds_DT <- as.data.frame(all_ds_DT)
  stopifnot(apply(all_ds_DT, 2, is.numeric))
  colnames(all_ds_DT) <- c(
    "nbrTADs",
    paste0("nbrTADsSignif_adjPvalThresh", all_pval_thresh),
    "nbrGenes",
    paste0("nbrGenesSignif_adjPvalThresh", all_pval_thresh),
    "ratioAUC_FCC",
    "ratioAUC_coexprDist")  
  stopifnot(apply(all_ds_DT, 2, is.numeric))
  all_ds_DT <- all_ds_DT[order(all_ds_DT$ratioAUC_FCC, decreasing=T),]
  all_ds_DT$ratioAUC_FCC<- round(all_ds_DT$ratioAUC_FCC, 4)
  all_ds_DT$ratioAUC_coexprDist <- round(all_ds_DT$ratioAUC_coexprDist, 4)
  all_ds_DT <- cbind(data.frame(dataset=rownames(all_ds_DT), stringsAsFactors = FALSE),
                     all_ds_DT)
  
  outFile <- file.path(outFold, "nbr_signif_pvalThresh_all_datasets.txt")
  write.table(all_ds_DT, file = outFile, sep="\t", quote=F, col.names=T, row.names=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

rownames(all_ds_DT) <- NULL
all_ds_DT <- all_ds_DT[order(all_ds_DT$ratioAUC_FCC, decreasing = TRUE),]
all_ds_DT$dataset <- factor(as.character(all_ds_DT$dataset), levels = as.character(all_ds_DT$dataset))


for(pval in all_pval_thresh) {
  all_ds_DT[, paste0("ratioTADsSignif_adjPvalThresh", pval)] <- round(all_ds_DT[, paste0("nbrTADsSignif_adjPvalThresh", pval)]/all_ds_DT$nbrTADs, 4)
  all_ds_DT[,paste0("ratioGenesSignif_adjPvalThresh",pval)] <- round(all_ds_DT[,paste0("nbrGenesSignif_adjPvalThresh", pval)]/all_ds_DT$nbrGenes, 4)
}


all_signifLevel <- c("all", "subset", "withTot")
all_signifType  <- c("Genes", "TADs")
all_yaxis <- c("ratio", "nbr")

signifLevel = "withTot"
signifType = "Genes"
yaxis = "nbr"

all_scatter_var_x <- colnames(all_ds_DT)[grepl("ratio.+adjPvalThresh.+", colnames(all_ds_DT))]
all_scatter_var_y <- c("ratioAUC_FCC", "ratioAUC_coexprDist")

my_plot_function <- function(var1, var2, DT, mylabels = NULL, withLab=F) {
  myx <- DT[,var1]
  myy <- DT[,var2]
  plot(x = myx,
       y = myy,
       xlab = paste0(var1),
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

myxvar = all_scatter_var_x[1]
myyvar = all_scatter_var_y[1]
for(myxvar in all_scatter_var_x) {
  for(myyvar in all_scatter_var_y){
    outFile <- file.path(outFold, paste0(myyvar, "_vs_", myxvar, "scatter.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    my_plot_function(var1=myxvar, var2=myyvar, DT=all_ds_DT, mylabels=all_ds_DT$dataset, withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

all_scatter_var_x <- colnames(all_ds_DT)[grepl("ratio.+adjPvalThresh.+", colnames(all_ds_DT))]
all_scatter_var_y <- all_scatter_var_x
for(myxvar in all_scatter_var_x) {
  for(myyvar in all_scatter_var_y){
    outFile <- file.path(outFold, paste0(myyvar, "_vs_", myxvar, "scatter.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    my_plot_function(var1=myxvar, var2=myyvar, DT=all_ds_DT, mylabels=all_ds_DT$dataset, withLab=T)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}



all_cols <- c("darkgreen", "darkorange", "darkorchid", "deeppink", "dodgerblue4", "firebrick", "gold", "slateblue", "darkgrey")
# all_cols <- c("dodgerblue4", "darkorange2")
stopifnot(length(all_cols) >= length(all_pval_thresh)+1)
  
for(signifLevel in all_signifLevel){
  for(signifType in all_signifType){
    for(yaxis in all_yaxis){
      
      if(signifLevel == "all") {
        pval_to_plot <- all_pval_thresh
      } else if(signifLevel == "subset") {
        pval_to_plot <- sub_pval_thresh
      }
      stopifnot(length(pval_to_plot) > 0)
      
      if(signifLevel == "withTot"){
        mycol <- all_cols[1:(length(pval_to_plot)+1)]
        pval_to_plot <- c(pval_to_plot,1)
      } else {
        mycol <- all_cols[1:length(pval_to_plot)]  
      }
      
      
      if(signifLevel == "withTot") {
        if(yaxis == "ratio") next
        signif_DT <- all_ds_DT[, grep(paste0(yaxis, signifType, "Signif_.+", "|dataset|nbr", signifType), colnames(all_ds_DT))]
        colnames(signif_DT)[colnames(signif_DT) == paste0("nbr", signifType)] <- paste0("nbr", signifType, "Signif_adjPvalThresh", 1)
      } else {
        signif_DT <- all_ds_DT[, grep(paste0(yaxis, signifType, "Signif_.+", "|dataset"), colnames(all_ds_DT))]
      }
      
      stopifnot(nrow(signif_DT) > 0)
      
      
      signif_DT_m <- melt(signif_DT, id="dataset")
      signif_DT_m$signifLevel <- gsub(paste0(".+_adjPvalThresh(.+)$"), "\\1", signif_DT_m$variable)
      # signif_DT_m$signifLevel[signif_DT_m$variable == paste0("nbr", signifType)] <- 1
      
      stopifnot(!is.na(as.numeric(as.character(signif_DT_m$signifLevel))))
      
      signif_DT_m <- signif_DT_m[as.numeric(as.character(signif_DT_m$signifLevel)) %in% pval_to_plot,]
      stopifnot(nrow(signif_DT_m) > 0)
      stopifnot(length(unique(signif_DT_m$signifLevel)) == length(unique(pval_to_plot)))
      signif_DT_m$signifLevel <- factor(signif_DT_m$signifLevel, 
                                        levels = as.character(sort(as.numeric(as.character(unique(signif_DT_m$signifLevel))), decreasing=T)))
      
      
      signifTypeTit <- ifelse(signifType != "TADs", tolower(signifType), signifType)
      
      plotTit <- paste0(toTitleCase(yaxis), " signif. ", signifTypeTit)
      mySub <- "(datasets ranked by ratioAUC FCC)"
      
      if(yaxis == "nbr") {
        signif_DT_m$value <- log10(signif_DT_m$value+1)
        ylab <- paste(yaxis, signifTypeTit, "signif. [log10(#+1)]")
      } else{
        ylab <- paste(yaxis, signifTypeTit, "signif.")
      }      
      p_signif <- ggplot(signif_DT_m, aes(x = dataset, y = value, fill = signifLevel)) +
        ggtitle(plotTit, subtitle = mySub)+
        geom_bar(stat="identity", position = "dodge")+
        scale_x_discrete(name="")+
        scale_y_continuous(name=ylab,
                           breaks = scales::pretty_breaks(n = 10))+
        scale_fill_manual(values = mycol) +
        labs(fill = "p-val. signif.\nthresh:")+
        theme( # Increase size of axis lines
          # top, right, bottom and left
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          plot.title = element_text(hjust = 0.5, face = "bold", size=16),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", size=10),
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
      if(SSHFS) p_signif
      
      outFile <- file.path(outFold, paste0(yaxis, "Signif_", signifType, "_", signifLevel, "PvalTresh", ".", plotType))
      ggsave(plot = p_signif, filename = outFile, height=myHeightGG, width = myWidthGG)
      cat(paste0("... written: ", outFile, "\n"))
      
    }
  }
}


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
