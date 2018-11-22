startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GOtargeted.R\n"))

# Rscript cmp_DE_TADs_genes_GOtargeted.R

# DA TADs genes vs. DE genes: C) how targeted are the results
# - number of different GO slims
# - entropy of the set of GO (kind of intra-set IC)

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ontologySimilarity, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  #GO_IC data
suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology, go, minimal_set
suppressPackageStartupMessages(library(GSEABase, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim

data(GO_IC) # GSEABase
data(go) # ontologyIndex

SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
source("analysis_utils.R")

topPercent <- 0.1
ontologyType <- "MF"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  if(args[1] %in% c("BP", "MF")){
    ontologyType <- args[1]
  }
}

outFold <- file.path("CMP_TADs_GENES_GOtargeted",ontologyType, topPercent)
# if(!SSHFS) system(paste0("mkdir -p ", outFold))
system(paste0("mkdir -p ", outFold))

# to retrieve FCC and coexprDist
dsFold <- file.path(setDir, 
                    "mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                    "OUTPUT_FOLDER")
stopifnot(file.exists(dsFold))

logFile <- file.path(outFold, "cmp_tads_genes_gotargeted_pvalselect_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

goSlimFile <- file.path(setDir, "/mnt/etemp/marie/GO_data/goslim_generic.obo")
goslimData <- get_ontology(goSlimFile) # ontologyIndex
GOslim <- getOBOCollection(goSlimFile)  # GSEABase

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... topPercent:\t", topPercent, "\n")
printAndLog(txt, logFile)
txt <- paste0("... ontologyType:\t", ontologyType, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:10]
# all_ds <- all_ds[1]
curr_ds=all_ds[1]
if(buildTable){
  # all_ds_DT <- foreach(curr_ds = all_ds, .combine='rbind') %dopar% {
    all_ds_DT <- foreach(curr_ds = all_ds, .combine='rbind') %do% {
      
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
    nTopTADs <- length(topTADs)
    
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
    nTopTADs_genes <- length(topTADs_genes)
    
    # take as many top genes as genes from topTADs  
    topGenes_manyAsTopTADs <- names(entrez_pval_sort[1:length(topTADs_genes)])
    
    nTopGenes_manyAsTopTADs <- length(topGenes_manyAsTopTADs)
    nTopGenes_manyAsPercent <- length(topGenes_manyAsTopPercent)
    
    stopifnot(length(topGenes_manyAsTopTADs) == length(topTADs_genes) )
    stopifnot(topTADs_genes %in% gene2tad_DT$entrezID)
    stopifnot(topTADs_genes %in% pipelineGenes)
    stopifnot(topGenes_manyAsTopPercent %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsTopPercent %in% pipelineGenes)
    stopifnot(topGenes_manyAsTopTADs %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsTopTADs %in% pipelineGenes)
    
    ##################
    ###### # of diff. go slim - topTADs_genes
    ##################
    # retrieve GO categories of the topTADs_genes

    topTADs_genes_GO_mapping <- bitr(topTADs_genes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db")

    topTADs_genes_GO_mapping <- topTADs_genes_GO_mapping[topTADs_genes_GO_mapping$ONTOLOGY == ontologyType,]
    topTADs_genes_GO <- unique(topTADs_genes_GO_mapping$GO)
    topTADs_genes_nGO <- length(topTADs_genes_GO)
    
    topTADs_genes_GOmin <- minimal_set(go, topTADs_genes_GO)
    topTADs_genes_nGOmin <- length(topTADs_genes_GOmin)
    
    
    # map the GO to their GO slim
    topTADs_genes_GO_collection <- GSEABase::GOCollection(topTADs_genes_GO)#GSEABase
    
    topTADs_genes_GO_slim <- goSlim(topTADs_genes_GO_collection, GOslim, ontologyType)
    
    topTADs_genes_GOslimFreq <-  setNames(topTADs_genes_GO_slim$Percent/100, rownames(topTADs_genes_GO_slim))
    stopifnot(abs(1-sum(topTADs_genes_GOslimFreq)) < 10^-6)
    txt <- paste0("... GOslim from the GO: ",  sum(topTADs_genes_GOslimFreq > 0), "/", length(topTADs_genes_GOslimFreq), "\n")
    printAndLog(txt, logFile)
    
    topTADs_genes_nGOslim <- length(topTADs_genes_GOslimFreq)
    
    # count the number of distinct go slim
    topTADs_genes_GOslimFreq <- topTADs_genes_GOslimFreq[topTADs_genes_GOslimFreq > 0]
    topTADs_genes_nGOslim <- length(topTADs_genes_GOslimFreq)
    
    # for each of the GO slim, compute their IC
    topTADs_genes_GOslim_IC <- -log2(topTADs_genes_GOslimFreq)
    stopifnot(!is.na(topTADs_genes_GOslim_IC))
    stopifnot(!is.infinite(topTADs_genes_GOslim_IC))
    topTADs_genes_GOslim_IC_mean <- mean(topTADs_genes_GOslim_IC)
    stopifnot(!is.na(topTADs_genes_GOslim_IC_mean))
    stopifnot(!is.infinite(topTADs_genes_GOslim_IC_mean))
    
    ##################
    ###### entropy GOs - topTADs_genes
    ##################
    # for all the GO categories from topTADs_genes, how often they pop up
    topTADs_genes_GOcount <- setNames(as.numeric(table(topTADs_genes_GO_mapping$GO)),
                                         names(table(topTADs_genes_GO_mapping$GO)))
    topTADs_genes_GOfreq <- topTADs_genes_GOcount/sum(topTADs_genes_GOcount) 
    # stopifnot(sum(topTADs_genes_GOfreq) == 1)
    topTADs_genes_GO_IC <- -log2(topTADs_genes_GOfreq)
    stopifnot(!is.na(topTADs_genes_GO_IC))
    stopifnot(!is.infinite(topTADs_genes_GO_IC))
    topTADs_genes_GO_IC_mean <- mean(topTADs_genes_GO_IC)
    stopifnot(!is.na(topTADs_genes_GO_IC_mean))
    stopifnot(!is.infinite(topTADs_genes_GO_IC_mean))
    
  ##################
  ###### # of diff. go slim - topGenes_manyAsTopPercent
  ##################
  # retrieve GO categories of the topGenes

  topGenes_manyAsTopPercent_GO_mapping <- bitr(topGenes_manyAsTopPercent, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db")

  topGenes_manyAsTopPercent_GO_mapping <- topGenes_manyAsTopPercent_GO_mapping[topGenes_manyAsTopPercent_GO_mapping$ONTOLOGY == ontologyType,]
  topGenes_manyAsTopPercent_GO <- unique(topGenes_manyAsTopPercent_GO_mapping$GO)
  topGenes_manyAsTopPercent_nGO <- length(topGenes_manyAsTopPercent_GO)
  
  topGenes_manyAsTopPercent_GOmin <- minimal_set(go, topGenes_manyAsTopPercent_GO)
  topGenes_manyAsTopPercent_nGOmin <- length(topGenes_manyAsTopPercent_GOmin)
  
  
  # map the GO to their GO slim
  topGenes_manyAsTopPercent_GO_collection <- GOCollection(topGenes_manyAsTopPercent_GO)#GSEABase
  topGenes_manyAsTopPercent_GO_slim <- goSlim(topGenes_manyAsTopPercent_GO_collection, GOslim, ontologyType)
  topGenes_manyAsTopPercent_GOslimFreq <-  setNames(topGenes_manyAsTopPercent_GO_slim$Percent/100, rownames(topGenes_manyAsTopPercent_GO_slim))
  stopifnot(abs(1-sum(topGenes_manyAsTopPercent_GOslimFreq)) < 10^-6)
  txt <- paste0("... GOslim from the GO: ",  sum(topGenes_manyAsTopPercent_GOslimFreq > 0), "/", length(topGenes_manyAsTopPercent_GOslimFreq), "\n")
  printAndLog(txt, logFile)
  
  topGenes_manyAsTopPercent_nGOslim <- length(topGenes_manyAsTopPercent_GOslimFreq)
  
  # count the number of distinct go slim
  topGenes_manyAsTopPercent_GOslimFreq <- topGenes_manyAsTopPercent_GOslimFreq[topGenes_manyAsTopPercent_GOslimFreq > 0]
  topGenes_manyAsTopPercent_nGOslim <- length(topGenes_manyAsTopPercent_GOslimFreq)
  
  # for each of the GO slim, compute their IC
  topGenes_manyAsTopPercent_GOslim_IC <- -log2(topGenes_manyAsTopPercent_GOslimFreq)
  stopifnot(!is.na(topGenes_manyAsTopPercent_GOslim_IC))
  stopifnot(!is.infinite(topGenes_manyAsTopPercent_GOslim_IC))
  topGenes_manyAsTopPercent_GOslim_IC_mean <- mean(topGenes_manyAsTopPercent_GOslim_IC)
  stopifnot(!is.na(topGenes_manyAsTopPercent_GOslim_IC_mean))
  stopifnot(!is.infinite(topGenes_manyAsTopPercent_GOslim_IC_mean))
  
  ##################
  ###### entropy GOs - topGenes_manyAsTopPercent
  ##################
  # for all the GO categories from topGenes_manyAsTopPercent, how often they pop up
  topGenes_manyAsTopPercent_GOcount <- setNames(as.numeric(table(topGenes_manyAsTopPercent_GO_mapping$GO)),
                                  names(table(topGenes_manyAsTopPercent_GO_mapping$GO)))
  topGenes_manyAsTopPercent_GOfreq <- topGenes_manyAsTopPercent_GOcount/sum(topGenes_manyAsTopPercent_GOcount) 
  # stopifnot(sum(topGenes_manyAsTopPercent_GOfreq) == 1)
  topGenes_manyAsTopPercent_GO_IC <- -log2(topGenes_manyAsTopPercent_GOfreq)
  stopifnot(!is.na(topGenes_manyAsTopPercent_GO_IC))
  stopifnot(!is.infinite(topGenes_manyAsTopPercent_GO_IC))
  topGenes_manyAsTopPercent_GO_IC_mean <- mean(topGenes_manyAsTopPercent_GO_IC)
  stopifnot(!is.na(topGenes_manyAsTopPercent_GO_IC_mean))
  stopifnot(!is.infinite(topGenes_manyAsTopPercent_GO_IC_mean))
  
  ##################
  ###### # of diff. go slim - topGenes_manyAsTopTADs
  ##################
  # retrieve GO categories of the topGenes

  topGenes_manyAsTopTADs_GO_mapping <- bitr(topGenes_manyAsTopTADs, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db")
  
  topGenes_manyAsTopTADs_GO_mapping <- topGenes_manyAsTopTADs_GO_mapping[topGenes_manyAsTopTADs_GO_mapping$ONTOLOGY == ontologyType,]
  topGenes_manyAsTopTADs_GO <- unique(topGenes_manyAsTopTADs_GO_mapping$GO)
  topGenes_manyAsTopTADs_nGO <- length(topGenes_manyAsTopTADs_GO)
  
  topGenes_manyAsTopTADs_GOmin <- minimal_set(go, topGenes_manyAsTopTADs_GO)
  topGenes_manyAsTopTADs_nGOmin <- length(topGenes_manyAsTopTADs_GOmin)
  
  
  # map the GO to their GO slim
  topGenes_manyAsTopTADs_GO_collection <- GOCollection(topGenes_manyAsTopTADs_GO)#GSEABase
  topGenes_manyAsTopTADs_GO_slim <- goSlim(topGenes_manyAsTopTADs_GO_collection, GOslim, ontologyType)
  topGenes_manyAsTopTADs_GOslimFreq <-  setNames(topGenes_manyAsTopTADs_GO_slim$Percent/100, rownames(topGenes_manyAsTopTADs_GO_slim))
  stopifnot(abs(1-sum(topGenes_manyAsTopTADs_GOslimFreq)) < 10^-6)
  txt <- paste0("... GOslim from the GO: ",  sum(topGenes_manyAsTopTADs_GOslimFreq > 0), "/", length(topGenes_manyAsTopTADs_GOslimFreq), "\n")
  printAndLog(txt, logFile)
  
  topGenes_manyAsTopTADs_nGOslim <- length(topGenes_manyAsTopTADs_GOslimFreq)
  
  # count the number of distinct go slim
  topGenes_manyAsTopTADs_GOslimFreq <- topGenes_manyAsTopTADs_GOslimFreq[topGenes_manyAsTopTADs_GOslimFreq > 0]
  topGenes_manyAsTopTADs_nGOslim <- length(topGenes_manyAsTopTADs_GOslimFreq)
  
  # for each of the GO slim, compute their IC
  topGenes_manyAsTopTADs_GOslim_IC <- -log2(topGenes_manyAsTopTADs_GOslimFreq)
  stopifnot(!is.na(topGenes_manyAsTopTADs_GOslim_IC))
  stopifnot(!is.infinite(topGenes_manyAsTopTADs_GOslim_IC))
  topGenes_manyAsTopTADs_GOslim_IC_mean <- mean(topGenes_manyAsTopTADs_GOslim_IC)
  stopifnot(!is.na(topGenes_manyAsTopTADs_GOslim_IC_mean))
  stopifnot(!is.infinite(topGenes_manyAsTopTADs_GOslim_IC_mean))
  
  ##################
  ###### entropy GOs - topGenes_manyAsTopTADs
  ##################
  # for all the GO categories from topGenes_manyAsTopTADs, how often they pop up
  topGenes_manyAsTopTADs_GOcount <- setNames(as.numeric(table(topGenes_manyAsTopTADs_GO_mapping$GO)),
                                            names(table(topGenes_manyAsTopTADs_GO_mapping$GO)))
  topGenes_manyAsTopTADs_GOfreq <- topGenes_manyAsTopTADs_GOcount/sum(topGenes_manyAsTopTADs_GOcount) 
  # stopifnot(sum(topGenes_manyAsTopTADs_GOfreq) == 1)
  topGenes_manyAsTopTADs_GO_IC <- -log2(topGenes_manyAsTopTADs_GOfreq)
  stopifnot(!is.na(topGenes_manyAsTopTADs_GO_IC))
  stopifnot(!is.infinite(topGenes_manyAsTopTADs_GO_IC))
  topGenes_manyAsTopTADs_GO_IC_mean <- mean(topGenes_manyAsTopTADs_GO_IC)
  stopifnot(!is.na(topGenes_manyAsTopTADs_GO_IC_mean))
  stopifnot(!is.infinite(topGenes_manyAsTopTADs_GO_IC_mean))
  
    data.frame(
      dataset=curr_ds,
      
      nTopTADs_genes = nTopTADs_genes,
      
      topTADs_genes_nGO=topTADs_genes_nGO,
      topTADs_genes_GO_IC_mean=topTADs_genes_GO_IC_mean,
      
      topTADs_genes_nGOmin=topTADs_genes_nGOmin,
      
      topTADs_genes_nGOslim=topTADs_genes_nGOslim,
      topTADs_genes_GOslim_IC_mean=topTADs_genes_GOslim_IC_mean,
      
      nTopGenes_manyAsPercent = nTopGenes_manyAsPercent,
      
      topGenes_manyAsTopPercent_nGO=topGenes_manyAsTopPercent_nGO,
      topGenes_manyAsTopPercent_GO_IC_mean=topGenes_manyAsTopPercent_GO_IC_mean,
      
      topGenes_manyAsTopPercent_nGOmin=topGenes_manyAsTopPercent_nGOmin,
      
      topGenes_manyAsTopPercent_nGOslim=topGenes_manyAsTopPercent_nGOslim,
      topGenes_manyAsTopPercent_GOslim_IC_mean=topGenes_manyAsTopPercent_GOslim_IC_mean,
      
      nTopGenes_manyAsTopTADs =nTopGenes_manyAsTopTADs,
      
      topGenes_manyAsTopTADs_nGO=topGenes_manyAsTopTADs_nGO,
      topGenes_manyAsTopTADs_GO_IC_mean=topGenes_manyAsTopTADs_GO_IC_mean,
      
      topGenes_manyAsTopTADs_nGOmin=topGenes_manyAsTopTADs_nGOmin,
      
      topGenes_manyAsTopTADs_nGOslim=topGenes_manyAsTopTADs_nGOslim,
      topGenes_manyAsTopTADs_GOslim_IC_mean=topGenes_manyAsTopTADs_GOslim_IC_mean,
      
      stringsAsFactors = FALSE
    )
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)

  all_ds_DT <- all_ds_DT[order(all_ds_DT$topTADs_genes_GO_IC_mean, 
                               all_ds_DT$topGenes_manyAsTopPercent_GO_IC_mean, 
                               all_ds_DT$topGenes_manyAsTopTADs_GO_IC_mean, 
                               decreasing=T),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}


all_ds_DT$topTADs_genes_nGO_ratio <- all_ds_DT$topTADs_genes_nGO/all_ds_DT$nTopTADs_genes
all_ds_DT$topTADs_genes_nGOslim_ratio <- all_ds_DT$topTADs_genes_nGOslim/all_ds_DT$nTopTADs_genes
all_ds_DT$topTADs_genes_nGOmin_ratio <- all_ds_DT$topTADs_genes_nGOmin/all_ds_DT$nTopTADs_genes

all_ds_DT$topGenes_manyAsTopPercent_nGO_ratio <- all_ds_DT$topGenes_manyAsTopPercent_nGO/all_ds_DT$nTopGenes_manyAsPercent
all_ds_DT$topGenes_manyAsTopPercent_nGOslim_ratio <- all_ds_DT$topGenes_manyAsTopPercent_nGOslim/all_ds_DT$nTopGenes_manyAsPercent
all_ds_DT$topGenes_manyAsTopPercent_nGOmin_ratio <- all_ds_DT$topGenes_manyAsTopPercent_nGOmin/all_ds_DT$nTopGenes_manyAsPercent

all_ds_DT$topGenes_manyAsTopTADs_nGO_ratio <- all_ds_DT$topGenes_manyAsTopTADs_nGO/all_ds_DT$nTopGenes_manyAsTopTADs
all_ds_DT$topGenes_manyAsTopTADs_nGOslim_ratio <- all_ds_DT$topGenes_manyAsTopTADs_nGOslim/all_ds_DT$nTopGenes_manyAsTopTADs
all_ds_DT$topGenes_manyAsTopTADs_nGOmin_ratio <- all_ds_DT$topGenes_manyAsTopTADs_nGOmin/all_ds_DT$nTopGenes_manyAsTopTADs

# stop("-- ok --\n")

# colnames(all_ds_DT)
# [1] "dataset"                                  "topTADs_genes_nGO"                       
# [3] "topTADs_genes_GO_IC_mean"                 "topTADs_genes_nGOslim"                   
# [5] "topTADs_genes_GOslim_IC_mean"             "topGenes_manyAsTopPercent_nGO"           
# [7] "topGenes_manyAsTopPercent_GO_IC_mean"     "topGenes_manyAsTopPercent_nGOslim"       
# [9] "topGenes_manyAsTopPercent_GOslim_IC_mean" "topGenes_manyAsTopTADs_nGO"              
# [11] "topGenes_manyAsTopTADs_GO_IC_mean"        "topGenes_manyAsTopTADs_nGOslim"          
# [13] "topGenes_manyAsTopTADs_GOslim_IC_mean"   

######################################################################## retrieve FCC and coexprDist

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

plotDT_m <- melt(all_ds_DT, id="dataset")

plotDT_m$dataset <- factor(plotDT_m$dataset, levels=dataset_order_FCC)

######################################################################## 


all_vars_toplot <- c("GO_IC_mean", "GOslim_IC_mean", "nGO", "nGOmin","nGOslim",  "nGO_ratio", "nGOmin_ratio","nGOslim_ratio")
ref_var <- "topTADs_genes"
all_cmp_var <- c("topGenes_manyAsTopPercent", "topGenes_manyAsTopTADs")


var_to_plot=all_vars_toplot[1]
cmp_var=all_cmp_var[1]

for(var_to_plot in all_vars_toplot){
  for(cmp_var in all_cmp_var){
    
    var1 <- paste0(ref_var,"_",var_to_plot)
    var2 <- paste0(cmp_var,"_",var_to_plot)
    
    curr_plotDT <- plotDT_m[plotDT_m$variable %in% c(var1,var2 ),]
    
    
    plotTit <- gsub("_", " ", var_to_plot)
    mySub <- paste0("ontology = ", ontologyType, 
                    " - topPercent = ", topPercent,
                    " - dataset order = FCC")
    
    
    curr_plotDT$variable <- gsub(paste0("_", var_to_plot), "", curr_plotDT$variable)
    
    
    p_var <- ggplot(curr_plotDT, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), c(ref_var, cmp_var)),
                        labels = setNames(c(ref_var, cmp_var), c(ref_var, cmp_var)))+
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
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_", ref_var, "_", cmp_var, "_barplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    curr_plotDT$variable <- as.character(curr_plotDT$variable)
    
    my_comparisons <- list(c(ref_var, cmp_var))
    p_box <- ggviolin(curr_plotDT, x = "variable", y = "value", 
             fill = "variable",
             palette = c("#00AFBB", "#FC4E07"),
             title = plotTit,
             xlab="",
             ylab = gsub("_", " ", var_to_plot),
             # sub=mySub,
             legend.title = "",
             add = "boxplot", add.params = list(fill = "white"))+
      stat_compare_means(comparisons = my_comparisons, 
                         # position = "bottomleft",
                         label = "p.signif")+ # Add significance levels
      stat_compare_means()+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_", ref_var, "_", cmp_var, "_boxplot.", plotType))
    ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}

all_vars_toplot <- c("GO_IC_mean", "GOslim_IC_mean", "nGO", "nGOmin","nGOslim", "nGO_ratio", "nGOmin_ratio","nGOslim_ratio")
all_ref_var <- c("aucFCC", "aucCoexprDist")
all_cmp_var <- c("topTADs_genes","topGenes_manyAsTopPercent", "topGenes_manyAsTopTADs")

var_to_plot <- all_vars_toplot[1]
ref_var <- all_ref_var[1]
cmp_var <- all_cmp_var[1]


for(var_to_plot in all_vars_toplot) {
  for(ref_var in all_ref_var){
    for(cmp_var in all_cmp_var) {
      
      var1 <- ref_var
      var2 <- paste0(cmp_var, "_", var_to_plot)
      
      stopifnot(var1 %in% colnames(all_ds_DT))
      stopifnot(var2 %in% colnames(all_ds_DT))
      
      myx <- all_ds_DT[,var2]
      myy <- all_ds_DT[,var1]
      outFile <- file.path(outFold, paste0(ref_var, "_", var2, "_scatterplot.", plotType))
      do.call(plotType, list(outFile, height = myHeight, width = myWidth))
      plot(x=myx,
           y=myy,
           main = paste0(gsub("_", " ", var_to_plot), " vs. ", ref_var),
           # sub = paste0(cmp_var),
           xlab=paste0(gsub("_", " ", var2)),
           ylab=var1,
           pch=16,cex=0.7)
      mtext(cmp_var, side=3)
      text(x=myx,y=myy,cex=0.7,labels=all_ds_DT[,"dataset"])
      add_curv_fit(x=myx,
                   y=myy, withR2 = F, lty=2)
      addCorr(x=myx,
              y=myy, legPos = "topleft", bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
    }
  }
}


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
