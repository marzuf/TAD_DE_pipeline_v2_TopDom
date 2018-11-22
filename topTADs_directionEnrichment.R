startTime <- Sys.time()
cat(paste0("> Rscript topTADs_directionEnrichment.R\n"))

# Rscript topTADs_directionEnrichment.R


options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- TRUE

source("analysis_utils.R")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 500, 10)


registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

outFold <- "TOP_TADS_ENRICH_DIR"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "topTADs_enrich_dir_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)



txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... topPercent:\t", topPercent, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:10]
# curr_ds = "TCGAcrc_msi_mss"
# curr_ds="GSE102073_stic_nostic"
# all_ds <- all_ds[1:3]

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
    
    ### RETRIEVE THE RATIO DOWN
    step8_fold <- file.path(dsFold, curr_ds, "8c_runAllDown")
    stopifnot(file.exists(step8_fold))
    ratioDownFile <- file.path(step8_fold, "all_obs_ratioDown.Rdata")
    stopifnot(file.exists(ratioDownFile))
    
    ratioDown <- eval(parse(text = load(ratioDownFile)))
    stopifnot(topTADs %in% names(ratioDown))
    ratioDown <- ratioDown[topTADs]
    
    ### RETRIEVE THE MEAN EXPRESSION
    settingF <- file.path(setDir,
                          "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput",
                          paste0("run_settings_", curr_ds, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)
    sample1_file <- file.path(setDir, sample1_file)
    stopifnot(file.exists(sample1_file))
    s1 <- eval(parse(text = load(sample1_file)))
    
    sample2_file <- file.path(setDir, sample2_file)
    stopifnot(file.exists(sample2_file))
    s2 <- eval(parse(text = load(sample2_file)))
    
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    
    rnaFile <- file.path(step0_fold, "rna_fpkmDT.Rdata")
    stopifnot(file.exists(rnaFile))
    rnaDT <- eval(parse(text = load(rnaFile)))
    
    # ratio of expression for the topTADs
    signExpr <- sapply(topTADs, function(x) {
      
      # retrieve the genes that are in the current TAD
      tad_genes <- as.character(gene2tad_DT$entrezID[as.character(gene2tad_DT$region) == as.character(x)])
      tad_genes <- tad_genes[tad_genes %in% pipelineGenes]  
      stopifnot(length(tad_genes) > 0)
      tad_genes_rownames <- names(pipelineGenes[pipelineGenes %in% tad_genes])
      stopifnot(length(tad_genes_rownames) > 0)
      stopifnot(!is.na(tad_genes_rownames))
      subDT1 <- as.numeric(unlist(rnaDT[tad_genes_rownames, s1]))
      stopifnot(length(subDT1) > 0)
      stopifnot(length(subDT1) == length(s1) * length(tad_genes_rownames))
      meanExpr1 <- mean(subDT1)
      stopifnot(!is.na(meanExpr1))
      subDT2 <- as.numeric(unlist(rnaDT[tad_genes_rownames, s2]))
      stopifnot(length(subDT2) > 0)
      stopifnot(length(subDT2) == length(s2) * length(tad_genes_rownames))
      meanExpr2 <- mean(subDT2)
      stopifnot(!is.na(meanExpr2))
      as.numeric(meanExpr1>meanExpr2)
      
    })
    # for the topTADs:
    # ratio of topTADs with ratioDown > 0.5
    # ratio of topTADs with signExpr == 1
    data.frame(
      dataset = curr_ds,
      nTopTADs = length(topTADs),
      signExprRatio = sum(signExpr == 1)/length(topTADs),
      ratioDownRatio= sum(ratioDown > 0.5)/length(topTADs),
      stringsAsFactors = FALSE
    )
  }
  outFile <- file.path(outFold, "all_ds_DT.Rdata" )
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata" )
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

sum( (all_ds_DT$signExprRatio > 0.5 & all_ds_DT$ratioDownRatio > 0.5) | (all_ds_DT$signExprRatio < 0.5 & all_ds_DT$ratioDownRatio < 0.5 ) )
# 60/66

outFile <- file.path(outFold, paste0("ratios_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  signExprRatio = all_ds_DT$signExprRatio,
  ratioDownRatio = all_ds_DT$ratioDownRatio
),
plotTit = "Ratio topTADs with\nTAD ratioDown > 0.5 or TAD meanExpr1 > meanExpr2",
my_xlab = "ratio"
  )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

