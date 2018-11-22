startTime <- Sys.time()
cat(paste0("> Rscript GOgraph_minimal_set_effect.R\n"))

# Rscript GOgraph_minimal_set_effect.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(AnnotationDbi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(topGO, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) #inducedGraph
suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# retrieved from ontologyIndex; needed for ontologyIndex::minimal_set()
data(go)

allGOs_AnnotationDBI <- keys(GO.db)
allGOs_ontologyIndex <- go$id

GOid_GOterm_DT_GOdb <- select(GO.db, keys=keys(GO.db), columns=c("GOID", "TERM"))
GOid_GOterm_DT_go <- data.frame(GOID = go$id, TERM=go$name, stringsAsFactors = F)
length(GOid_GOterm_DT_GOdb$TERM)
#44541
length(GOid_GOterm_DT_go$TERM)
# 45638
length(intersect(GOid_GOterm_DT_GOdb$TERM,GOid_GOterm_DT_go$TERM))
# 43442
intersectIDs <- intersect(GOid_GOterm_DT_GOdb$GOID,GOid_GOterm_DT_go$GOID)
length(intersectIDs)
# 43569
stopifnot(!any(duplicated(GOid_GOterm_DT_go$GOID)))
stopifnot(!any(duplicated(GOid_GOterm_DT_GOdb$GOID)))
# stopifnot(!any(duplicated(GOid_GOterm_DT_go$TERM)))
# stopifnot(!any(duplicated(GOid_GOterm_DT_GOdb$TERM)))

# they don't use the same terms (one uses abreviations in the terms, not the other)

GOid_GOterm_DT_cP <- clusterProfiler:::get_GO2TERM_table()
colnames(GOid_GOterm_DT_cP) <- c("GOID", "TERM")
stopifnot(!any(duplicated(GOid_GOterm_DT_cP$GOID)))
stopifnot(!any(duplicated(GOid_GOterm_DT_cP$TERM)))
GOterm_GOid_cP <- setNames(GOid_GOterm_DT_cP$GOID, GOid_GOterm_DT_cP$TERM)


length(allGOs_AnnotationDBI)
# 44541
length(allGOs_ontologyIndex)
# 45638
length(intersect(allGOs_ontologyIndex, allGOs_AnnotationDBI))
# 43569

SSHFS <- FALSE
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
source("GO_graph_utils.R")

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
topPercent <- 0.10
pvalSelect <- 0.05
pvalSelectGO <- 0.05

curr_ds <- "TCGAcrc_msi_mss"

enricher_ontologyType <- "BP"

outFold <- file.path("GOgraph_minimal_set_effect",  enricher_ontologyType, paste0(topPercent, "_", pvalSelect, "_", pvalSelectGO), curr_ds)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)


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

GO_g <- makeGOGraph(ont = tolower(enricher_ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)

enricher_results_sortGOby <- "p.adjust"


# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))


mySub <- paste0("(select top percent = ", topPercent , "; signif. GO = ",pvalSelectGO, ")") 

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
txt <- paste0("... enricher ontologyType:\t", enricher_ontologyType, "\n")
printAndLog(txt, logFile)
txt <- paste0("... topPercent select genes:\t", topPercent, "\n")
printAndLog(txt, logFile)
txt <- paste0("... pval thresh select genes:\t", pvalSelect, "\n")
printAndLog(txt, logFile)
txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:5]

c5_msigdb <- read.gmt(gmtFile)

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
topGenes_manyAsPercent <- names(entrez_pval_sort[1:nTopGenes])
nTopGenes_manyAsPercent <- length(topGenes_manyAsPercent)

stopifnot(nTopGenes > 0)
stopifnot(topTADs %in% gene2tad_DT$region)

topTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% topTADs]  
topTADs_genes <- topTADs_genes[topTADs_genes %in% pipelineGenes]
nTopTADs_genes <- length(topTADs_genes)
stopifnot(nTopTADs_genes > 0)
stopifnot(nTopTADs_genes <= length(entrez_pval_sort) )

# take as many top genes as genes from topTADs  
topGenes_manyAsTopTADs <- names(entrez_pval_sort[1:length(topTADs_genes)])
nTopGenes_manyAsTopTADs <- length(topGenes_manyAsTopTADs)

stopifnot(length(topGenes_manyAsTopTADs) == length(topTADs_genes) )
stopifnot(topTADs_genes %in% gene2tad_DT$entrezID)
stopifnot(topTADs_genes %in% pipelineGenes)
stopifnot(topGenes_manyAsPercent %in% gene2tad_DT$entrezID)
stopifnot(topGenes_manyAsPercent %in% pipelineGenes)
stopifnot(topGenes_manyAsTopTADs %in% gene2tad_DT$entrezID)
stopifnot(topGenes_manyAsTopTADs %in% pipelineGenes)

selectTADs <- names(tad_pval[tad_pval <= pvalSelect])
nSelectTADs <- length(selectTADs)


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
  nSelectTADs_genes <- length(selectTADs_genes)
  
} else {
  nSelectTADs_genes <- 0
}

if(nSelectGenes > 0){
  stopifnot(selectGenes %in% gene2tad_DT$entrezID)
  stopifnot(selectGenes %in% pipelineGenes)
}



############################################################################################
############################################################################################ selectTADs_genes
############################################################################################

#***** 1) selectTADs_genes

selectTADs_genes_enrich <- enricher(gene = selectTADs_genes, 
                                  TERM2GENE=c5_msigdb,
                                  pvalueCutoff = enricher_pvalueCutoff, 
                                  pAdjustMethod = enricher_pAdjustMethod, 
                                  minGSSize = enricher_minGSSize, 
                                  maxGSSize = enricher_maxGSSize, 
                                  qvalueCutoff =enricher_qvalueCutoff)

selectTADs_genes_resultDT <- selectTADs_genes_enrich@result
selectTADs_genes_resultDT <- selectTADs_genes_resultDT[order(selectTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

selectTADs_genes_signifGOterm <- as.character(selectTADs_genes_resultDT$ID[selectTADs_genes_resultDT$p.adjust <= pvalSelectGO])
nSelectTADs_genes_signifGOterm <- length(selectTADs_genes_signifGOterm)

selectTADs_genes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", selectTADs_genes_signifGOterm)))

length(selectTADs_genes_signifGOterm_tmp)
# 58
sum(selectTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
# 46
sum(selectTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
# 46
sum(selectTADs_genes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
# 46

selectTADs_genes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% selectTADs_genes_signifGOterm_tmp]
nSelectTADs_genes_signifGOid <- length(selectTADs_genes_signifGOid)


###### GO analysis

selectTADs_genes_signifGOmin <- minimal_set(go, selectTADs_genes_signifGOid)
nSelectTADs_genes_signifGOmin <- length(selectTADs_genes_signifGOmin)
selectTADs_genes_gNEL <- inducedGraph(dag = GO_g, startNodes = selectTADs_genes_signifGOmin[selectTADs_genes_signifGOmin %in% nodes(GO_g)])

# !!! need to reverse the graph !!!
selectTADs_genes_gNEL_rev <- reverseArch(selectTADs_genes_gNEL) # topGO
selectTADs_genes_gNEL_rev_noAll <- removeNode("all", selectTADs_genes_gNEL_rev) # not needed ? -> do not want this root
selectTADs_genes_g <- igraph.from.graphNEL(graphNEL = selectTADs_genes_gNEL_rev_noAll)

plot(selectTADs_genes_gNEL)

###### GO analysis - gene subset

addAncestor <- TRUE

for(i in 1:2){
  selectTADs_genes_signifGOsub <- selectTADs_genes_signifGOid[1:2]
  # add some of the parents to test the minimal set
  
  if(addAncestor) {
    GO_ANCESTOR <- as.list(eval(parse(text = paste0("GO", enricher_ontologyType, "ANCESTOR"))))
    someAncestors <- unique(unlist(GO_ANCESTOR[selectTADs_genes_signifGOsub]))
    someAncestors <- someAncestors[someAncestors != "all"]
    someAncestors <- sample(someAncestors, size = ceiling(length(someAncestors)/2))
    selectTADs_genes_signifGOsub <- c(selectTADs_genes_signifGOsub, someAncestors)
    outFile1 <- file.path(outFold, paste0("selectTADs_genes_beforeMin_withAncestors.", plotType))
    outFile2 <- file.path(outFold, paste0("selectTADs_genes_afterMin_withAncestors.", plotType))
  } else{
    outFile1 <- file.path(outFold, paste0("selectTADs_genes_beforeMin.", plotType))  
    outFile2 <- file.path(outFold, paste0("selectTADs_genes_afterMin.", plotType))  
  }
  
  GO_to_plot <- selectTADs_genes_signifGOsub[selectTADs_genes_signifGOsub %in% nodes(GO_g)]
  selectTADs_genes_gNELsub <- inducedGraph(dag = GO_g, startNodes = GO_to_plot)
  # !!! need to reverse the graph !!!
  selectTADs_genes_gNELsub_rev <- reverseArch(selectTADs_genes_gNELsub) # topGO
  selectTADs_genes_gNELsub_rev_noAll <- removeNode("all", selectTADs_genes_gNELsub_rev) # not needed ? -> do not want this root
  selectTADs_genes_gSub <- igraph.from.graphNEL(graphNEL = selectTADs_genes_gNELsub_rev_noAll)
  
  
  nAttrsSub<-list()
  nAttrsSub$color <- setNames(rep("black", length(nodes(selectTADs_genes_gNELsub_rev_noAll))), nodes(selectTADs_genes_gNELsub_rev_noAll))
  nAttrsSub$color[GO_to_plot] <- "red"
  nAttrsSub$fillcolor <- setNames(rep("transparent", length(nodes(selectTADs_genes_gNELsub_rev_noAll))), nodes(selectTADs_genes_gNELsub_rev_noAll))
  nAttrsSub$fillcolor[GO_to_plot] <- "orange"
  
  do.call(plotType, list(outFile1, height=myHeight, width=myWidth))
  plot(selectTADs_genes_gNELsub_rev_noAll, nodeAttrs = nAttrsSub)
  mtext(paste0("(before minimal set)\nGO to draw: ", paste0(GO_to_plot, collapse=", ")), side = 3, line=-2)
  foo <- dev.off()
  cat(paste0("... written: ", outFile1, "\n"))
  
  
  
  ###### GO analysis - gene subset - with minimal_set
  selectTADs_genes_signifGOsubMin <- minimal_set(go, selectTADs_genes_signifGOsub)
  GOmin_to_plot <- selectTADs_genes_signifGOsubMin[selectTADs_genes_signifGOsubMin %in% nodes(GO_g)]
  selectTADs_genes_gNELsubMin <- inducedGraph(dag = GO_g, startNodes = GOmin_to_plot)
  # !!! need to reverse the graph !!!
  selectTADs_genes_gNELsubMin_rev <- reverseArch(selectTADs_genes_gNELsubMin) # topGO
  selectTADs_genes_gNELsubMin_rev_noAll <- removeNode("all", selectTADs_genes_gNELsubMin_rev) # not needed ? -> do not want this root
  selectTADs_genes_gsubMin <- igraph.from.graphNEL(graphNEL = selectTADs_genes_gNELsubMin_rev_noAll)
  
  nAttrsSubMin<-list()
  nAttrsSubMin$color <- setNames(rep("blue", length(nodes(selectTADs_genes_gNELsubMin_rev_noAll))), nodes(selectTADs_genes_gNELsubMin_rev_noAll))
  nAttrsSubMin$color[GOmin_to_plot] <- "red"
  nAttrsSub$fillcolor <- setNames(rep("transparent", length(nodes(selectTADs_genes_gNELsub_rev_noAll))), nodes(selectTADs_genes_gNELsub_rev_noAll))
  nAttrsSub$fillcolor[GOmin_to_plot] <- "orange"
  
  do.call(plotType, list(outFile2, height=myHeight, width=myWidth))
  plot(selectTADs_genes_gNELsubMin_rev_noAll, nodeAttrs = nAttrsSubMin)
  mtext(paste0("(after minimal_set)\nGO to draw: ", paste0(GOmin_to_plot, collapse=", ")), side = 3, line=-2)
  foo <- dev.off()
  cat(paste0("... written: ", outFile2, "\n"))
  
  addAncestor <- !addAncestor
}








######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

