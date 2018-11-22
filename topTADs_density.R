startTime <- Sys.time()
cat(paste0("> Rscript topTADs_density.R\n"))

# Rscript topTADs_density.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

source("analysis_utils.R")

SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 500, 10)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

nTopTADs <- 20
# nTopTADs <- "all"
nToPlot <- 5

outFold <- file.path("~/Desktop/TOP_TADs_DENSITY",  nTopTADs)

# outFold <- file.path("TOP_TADs_DENSITY",  nTopTADs)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "top_tads_density.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

mySub <- paste0("(nTopTADs = ", nTopTADs , ")") 

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... nTopTADs:\t", nTopTADs, "\n")
printAndLog(txt, logFile)
txt <- paste0("... nToPlot:\t", nToPlot, "\n")
printAndLog(txt, logFile)


if(nTopTADs!="all") stopifnot(nTopTADs > 0)


# all_ds <- all_ds[1:2]


#all_ds <- all_ds[1:3]
curr_ds="TCGAcrc_msi_mss"
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
    tad_pval <- sort(tad_pval, decreasing=FALSE)
    
    
    if(nTopTADs == "all") {
      topTADs <- names(tad_pval)
    } else {
      topTADs <- names(tad_pval[1:nTopTADs])
    }
    
    
    topTADs_g2t_DT <- gene2tad_DT[gene2tad_DT$region  %in% topTADs,]
    stopifnot(nrow(topTADs_g2t_DT) > 0)
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(length(pipelineGenes) > 0)
    stopifnot(pipelineGenes %in% gene2tad_DT$entrezID)
    
    topTADs_g2t_DT <- topTADs_g2t_DT[as.character(topTADs_g2t_DT$entrezID) %in% pipelineGenes,]
    
    stopifnot(topTADs %in% topTADs_g2t_DT$region)

    topTADs_nbrGenes <- setNames(as.numeric(table(topTADs_g2t_DT$region)), names(table(topTADs_g2t_DT$region)))
    stopifnot(topTADs %in% names(topTADs_nbrGenes))
    
    data.frame(
      dataset = curr_ds, 
      region = topTADs,
      nbrGenes = topTADs_nbrGenes[topTADs],
      stringsAsFactors=FALSE
    )

    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

# stop("-- ok -- \n")

all_ds_aucFCC <- foreach(curr_ds = unique(all_ds_DT$dataset), .combine='c') %dopar% {
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
names(all_ds_aucFCC) <- unique(all_ds_DT$dataset)
all_ds_aucFCC <- sort(all_ds_aucFCC, decreasing = TRUE)

all_ds_aucCoexprDist <- foreach(curr_ds = unique(all_ds_DT$dataset), .combine='c') %dopar% {
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
names(all_ds_aucCoexprDist) <- unique(all_ds_DT$dataset)
all_ds_aucCoexprDist <- sort(all_ds_aucCoexprDist, decreasing = TRUE)

splittedNbrGenes <- split(all_ds_DT$nbrGenes, unique(all_ds_DT$dataset))

############################################################################################################
############################################################################################################ DENSITY PLOTS
############################################################################################################

topDs_FCC <- names(all_ds_aucFCC)[1:nToPlot]
botDs_FCC <- names(all_ds_aucFCC)[(length(all_ds_aucFCC) - nToPlot + 1):length(all_ds_aucFCC)]
stopifnot(length(topDs_FCC) == length(botDs_FCC) )

outFile <- file.path(outFold, paste0("multidens_topBotDs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(splittedNbrGenes[c(topDs_FCC, botDs_FCC)],
               plotTit="nbrGenes/TAD - topFCC and botFCC datasets", legTxt=NULL, legPos="topright",
               my_xlab="nbrGenes")
mtext(mySub, side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n" ))
  
outFile <- file.path(outFold, paste0("multidens_topDs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(splittedNbrGenes[c(topDs_FCC)],
               plotTit="nbrGenes/TAD - topFCC datasets", legTxt=NULL, legPos="topright",
               my_xlab="nbrGenes")
mtext(mySub, side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n" ))
  
outFile <- file.path(outFold, paste0("multidens_botDs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(splittedNbrGenes[c(botDs_FCC)],
               plotTit="nbrGenes/TAD - botFCC datasets", legTxt=NULL, legPos="topright",
               my_xlab="nbrGenes")
mtext(mySub, side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n" ))

cat("AAA\n")
############################################################################################################
############################################################################################################ MEAN/MEDIAN/QUANTILE ~ AUC 
############################################################################################################

# dt1 <- aggregate(nbrGenes ~ dataset, data=all_ds_DT, FUN=min)
# colnames(dt1)[2] = "min"
# dt2 <- aggregate(nbrGenes ~ dataset, data=all_ds_DT, FUN=mean)
# colnames(dt2)[2] = "mean"
# dt3 <- aggregate(nbrGenes ~ dataset, data=all_ds_DT, FUN=median)
# colnames(dt3)[2] = "median"
# dt4 <- aggregate(nbrGenes ~ dataset, data=all_ds_DT, FUN=max)
# colnames(dt4)[2] = "max"
# dt5 <- aggregate(nbrGenes ~ dataset, data=all_ds_DT, FUN=length)
# colnames(dt5)[2] = "nbr"
# 
# summaryDT <- merge(dt5, merge(dt4, merge(dt3, merge(dt1, dt2, by="dataset"),by ="dataset"), by="dataset"), by ="dataset")

# summaryDT <- by(all_ds_DT, all_ds_DT$dataset, function(x) {
#   # return(1)
#   # save(x, file="x.Rdata")
#   # tmp <- setNames(as.numeric(summary(x$nbrGenes)), as.character(names(summary(x$nbrGenes))))
#   # tmp <- c(tmp, Tot=sum(x$nbrGenes))
#   tmp <- c(a=1,b=2,c=3)
#   # tmp <- c(as.numeric(summary(x$nbrGenes)), sum(x$nbrGenes)) 
#   # names(tmp) <- c(as.character(names(summary(x$nbrGenes))), "tot")
#   # # as.character(names(summary(x$nbrGenes))))
#   # # tmp <- c(tmp, Tot=)
#   # names(tmp) <- gsub(" ", "", gsub("\\.", "", names(tmp)))
#   return(tmp)
# })
# stop("--ok\n")
# 
summaryDT <- as.data.frame(do.call(rbind, by(all_ds_DT, all_ds_DT$dataset, function(x) {
  tmp <- setNames(as.numeric(summary(x$nbrGenes)), as.character(names(summary(x$nbrGenes))))
  tmp <- c(tmp, Tot=sum(x$nbrGenes))#
  # tmp <- c(as.numeric(summary(x$nbrGenes)), sum(x$nbrGenes))
  # # as.character(names(summary(x$nbrGenes))))
  # # tmp <- c(tmp, Tot=)
  # names(tmp) <- c(gsub(" ", "", gsub("\\.", "", names(tmp))), "Tot")
  tmp
})))


summaryDT$dataset <- rownames(summaryDT)
summaryDT$dataset <- rownames(summaryDT)


summaryDT$aucFCC <- all_ds_aucFCC[summaryDT$dataset]
summaryDT$aucCoexprDist <- all_ds_aucCoexprDist[summaryDT$dataset]

cat("DDD\n")

mysub <- paste0("(nTopTADs = ", nTopTADs , ")") 

all_y_vars <- colnames(summaryDT)
all_x_vars <- c("aucFCC", "aucCoexprDist")
all_y_vars <- all_y_vars[all_y_vars!="dataset" & !all_y_vars %in%all_x_vars]

yvar=all_y_vars[1]
xvar=all_x_vars[1]
yvar="mean"
cat("CCC\n")


for(yvar in all_y_vars){
  for(xvar in all_x_vars){
    
    myylab <- paste0("# genes - ", yvar)
    myxlab <- paste0(xvar)
    
    mymain <- paste0(yvar, " # genes vs. ", xvar)
    
    myx <- summaryDT[,xvar]
    myy <- summaryDT[,yvar]
    
    outFile <- file.path(outFold, paste0(xvar, "_", gsub(" ", "", gsub("\\.", "", yvar)), ".", plotType))
    # outFile="foo.svg"
    # plot(NULL, xlim=c(0,1), ylim=c(0,1))
    # foo <- dev.off()
    
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x=myx,
         y=myy,
         main = mymain,
         xlab =  myxlab,
         ylab = myylab,
         pch=16,cex=0.7)
    # foo <- dev.off()
    add_curv_fit(x=myx,
                 y=myy, withR2 = F, lty=2)
    addCorr(x=myx,
            y=myy, legPos = "topleft", bty="n")
    text(x = myx, y = myy, labels = summaryDT[,"dataset"], cex=0.7)
    mtext(mysub,side=3)
    foo <- dev.off()
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

