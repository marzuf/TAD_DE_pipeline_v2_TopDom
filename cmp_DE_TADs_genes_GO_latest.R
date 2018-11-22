startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO_latest.R\n"))

# Rscript cmp_DE_TADs_genes_GO_latest.R

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

SSHFS <- F
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

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

pvalSelectGO <- 0.05

outFold <- "CMP_TADs_GENES_GO_latest"
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



# settings for enrichment analysis: clusterProfiler::enricher
# default values enricher function: (universe as background !)
# pvalueCutoff = 0.05 
# pAdjustMethod = "BH"
# minGSSize = 10
# maxGSSize = 500
# qvalueCutoff = 0.2 
enricher_ontologyType <- "BP"
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
txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:5]

c5_msigdb <- read.gmt(gmtFile)

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
  
    
    ############################################################################################
    ############################################################################################ enricher - topTADs
    ############################################################################################
    
  # run GO analysis for these 3 sets of entrez genes: topTADs_genes, topGenes_manyAsPercent
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
  
  # topTADs_genes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topTADs_genes_resultDT)))
  # if(topTADs_genes_nTop > 0) {
  #   topTADs_genes_topGo  <- as.character(topTADs_genes_resultDT$ID[1:topTADs_genes_nTop])
  # } else {
  #   topTADs_genes_topGo  <- character(0)  
  # }
  # stopifnot(length(topTADs_genes_topGo) == topTADs_genes_nTop)
  
  topTADs_genes_signifGOterm <- as.character(topTADs_genes_resultDT$ID[topTADs_genes_resultDT$p.adjust <= pvalSelectGO])
  nTopTADs_genes_signifGOterm <- length(topTADs_genes_signifGOterm)
  
  topTADs_genes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topTADs_genes_signifGOterm)))
  
  length(topTADs_genes_signifGOterm_tmp)
  # 58
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 46
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 46
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 46
  
  topTADs_genes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topTADs_genes_signifGOterm_tmp]
  nTopTADs_genes_signifGOid <- length(topTADs_genes_signifGOid)
  
  # if(nTopTADs_genes_signifGOterm > 0)
  #   stopifnot(nTopTADs_genes_signifGOid > 0) # not always TRUE (task 3)
  

  ############################################################################################
  ############################################################################################ enricher - manyAsPercent
  ############################################################################################
  
  
#***** 2) topGenes_manyAsPercent
  topGenes_manyAsPercent_enrich <- enricher(gene = topGenes_manyAsPercent, 
                                 TERM2GENE=c5_msigdb,
                                 pvalueCutoff = enricher_pvalueCutoff, 
                                 pAdjustMethod = enricher_pAdjustMethod, 
                                 minGSSize = enricher_minGSSize, 
                                 maxGSSize = enricher_maxGSSize, 
                                 qvalueCutoff =enricher_qvalueCutoff)
  
  topGenes_manyAsPercent_resultDT <- topGenes_manyAsPercent_enrich@result
  topGenes_manyAsPercent_resultDT <- topGenes_manyAsPercent_resultDT[order(topGenes_manyAsPercent_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  
  # topGenes_manyAsPercent_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsPercent_resultDT)))
  # if(topGenes_manyAsPercent_nTop > 0) {
  #   topGenes_manyAsPercent_topGo  <- as.character(topGenes_manyAsPercent_resultDT$ID[1:topGenes_manyAsPercent_nTop])
  # } else {
  #   topGenes_manyAsPercent_topGo  <- character(0)  
  # }
  # stopifnot(length(topGenes_manyAsPercent_topGo) == topGenes_manyAsPercent_nTop)
  
  topGenes_manyAsPercent_signifGOterm <- as.character(topGenes_manyAsPercent_resultDT$ID[topGenes_manyAsPercent_resultDT$p.adjust <= pvalSelectGO])
  nTopGenes_manyAsPercent_signifGOterm <- length(topGenes_manyAsPercent_signifGOterm)
  
  topGenes_manyAsPercent_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topGenes_manyAsPercent_signifGOterm)))
  
  length(topGenes_manyAsPercent_signifGOterm_tmp)
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 5
  
  topGenes_manyAsPercent_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topGenes_manyAsPercent_signifGOterm_tmp]
  nTopGenes_manyAsPercent_signifGOid <- length(topGenes_manyAsPercent_signifGOid)
  
  # if(nTopGenes_manyAsPercent_signifGOterm > 0)
  #   stopifnot(nTopGenes_manyAsPercent_signifGOid > 0) # not true task 22
  
  ############################################################################################
  ############################################################################################ enricher - manyAsTopTADs
  ############################################################################################
  
  
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
  
  # topGenes_manyAsTopTADs_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsTopTADs_resultDT)))
  # if(topGenes_manyAsTopTADs_nTop > 0) {
  #   topGenes_manyAsTopTADs_topGo  <- as.character(topGenes_manyAsTopTADs_resultDT$ID[1:topGenes_manyAsTopTADs_nTop])
  # } else {
  #   topGenes_manyAsTopTADs_topGo  <- character(0)  
  # }
  # stopifnot(length(topGenes_manyAsTopTADs_topGo) == topGenes_manyAsTopTADs_nTop)
  
  topGenes_manyAsTopTADs_signifGOterm <- as.character(topGenes_manyAsTopTADs_resultDT$ID[topGenes_manyAsTopTADs_resultDT$p.adjust <= pvalSelectGO])
  nTopGenes_manyAsTopTADs_signifGOterm <- length(topGenes_manyAsTopTADs_signifGOterm)
  
  topGenes_manyAsTopTADs_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topGenes_manyAsTopTADs_signifGOterm)))
  
  length(topGenes_manyAsTopTADs_signifGOterm_tmp)
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 5
  
  topGenes_manyAsTopTADs_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topGenes_manyAsTopTADs_signifGOterm_tmp]
  nTopGenes_manyAsTopTADs_signifGOid <- length(topGenes_manyAsTopTADs_signifGOid)
  
  # if(nTopGenes_manyAsPercent_signifGOterm > 0)
  #   stopifnot(nTopGenes_manyAsPercent_signifGOid > 0) # not true task 22
  
  
  
  ############################################################################################
  ############################################################################################ GRAPH - topTADs
  ############################################################################################
  

  ###### GO analysis
  
  if(!is.na(nTopTADs_genes_signifGOid) & nTopTADs_genes_signifGOid > 0) {
    
    topTADs_genes_signifGOmin <- minimal_set(go, topTADs_genes_signifGOid)
    nTopTADs_genes_signifGOmin <- length(topTADs_genes_signifGOmin)
    topTADs_genes_gNEL <- inducedGraph(dag = GO_g, startNodes = topTADs_genes_signifGOmin[topTADs_genes_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topTADs_genes_gNEL_rev <- reverseArch(topTADs_genes_gNEL) # topGO
    topTADs_genes_gNEL_rev_noAll <- removeNode("all", topTADs_genes_gNEL_rev) # not needed ? -> do not want this root
    topTADs_genes_g <- igraph.from.graphNEL(graphNEL = topTADs_genes_gNEL_rev_noAll)
    
    stopifnot(is.connected(topTADs_genes_g))
    stopifnot(is.directed(topTADs_genes_g))
    stopifnot(is.dag(topTADs_genes_g))
    
    topTADs_genes_g_nNodes <- length(V(topTADs_genes_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topTADs_genes_g_rootIdx <- which(sapply(sapply(V(topTADs_genes_g), 
    #                                 function(x) neighbors(topTADs_genes_g, x, mode="in")), length) == 0)
    # stopifnot(length(topTADs_genes_g_rootIdx) == 1)
    # topTADs_genes_g_rootGO <- names(V(topTADs_genes_g)[topTADs_genes_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topTADs_genes_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topTADs_genes_g_density <- edge_density(topTADs_genes_g)
    
    # edge connectivity (adhesion)
    # topTADs_genes_g_undir <- as.undirected(topTADs_genes_g)
    topTADs_genes_g_adhesion <- edge_connectivity(topTADs_genes_g)
    
    # vertex connectivity (cohesion)
    topTADs_genes_g_cohesion <- vertex_connectivity(topTADs_genes_g)
    # diameter
    topTADs_genes_g_diameter <- diameter(topTADs_genes_g)
    
    # mean geodesic distance
    topTADs_genes_g_meanDist <- mean_distance(topTADs_genes_g, directed = FALSE) # directed TRUE or FALSE ???
    topTADs_genes_g_meanDistDir <- mean_distance(topTADs_genes_g, directed = TRUE) # directed TRUE or FALSE ???
  
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topTADs_genes_g_meanEccentricity <- mean(eccentricity(topTADs_genes_g))

   # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topTADs_genes_g_meanBetweenness <- mean(betweenness(topTADs_genes_g))
    
  } else {
    topTADs_genes_signifGOmin <- NA
    nTopTADs_genes_signifGOmin <- NA
    topTADs_genes_g_nNodes <- NA
    topTADs_genes_g_density <- NA
    topTADs_genes_g_adhesion <- NA
    topTADs_genes_g_cohesion <- NA
    topTADs_genes_g_diameter <- NA
    topTADs_genes_g_meanDist <- NA
    topTADs_genes_g_meanDistDir <- NA
    topTADs_genes_g_meanEccentricity <- NA
    topTADs_genes_g_meanBetweenness <- NA
  }
    
  ############################################################################################
  ############################################################################################ GRAPH - manyAsPercent
  ############################################################################################
      
  if(!is.na(nTopGenes_manyAsPercent_signifGOid) & nTopGenes_manyAsPercent_signifGOid > 0) {
    
    topGenes_manyAsPercent_signifGOmin <- minimal_set(go, topGenes_manyAsPercent_signifGOid)
    nTopGenes_manyAsPercent_signifGOmin <- length(topGenes_manyAsPercent_signifGOmin)
    topGenes_manyAsPercent_gNEL <- inducedGraph(dag = GO_g, startNodes = topGenes_manyAsPercent_signifGOmin[topGenes_manyAsPercent_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topGenes_manyAsPercent_gNEL_rev <- reverseArch(topGenes_manyAsPercent_gNEL) # topGO
    topGenes_manyAsPercent_gNEL_rev_noAll <- removeNode("all", topGenes_manyAsPercent_gNEL_rev) # not needed ? -> do not want this root
    topGenes_manyAsPercent_g <- igraph.from.graphNEL(graphNEL = topGenes_manyAsPercent_gNEL_rev_noAll)
    
    stopifnot(is.connected(topGenes_manyAsPercent_g))
    stopifnot(is.directed(topGenes_manyAsPercent_g))
    stopifnot(is.dag(topGenes_manyAsPercent_g))
    
    topGenes_manyAsPercent_g_nNodes <- length(V(topGenes_manyAsPercent_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topGenes_manyAsPercent_g_rootIdx <- which(sapply(sapply(V(topGenes_manyAsPercent_g), 
    #                                 function(x) neighbors(topGenes_manyAsPercent_g, x, mode="in")), length) == 0)
    # stopifnot(length(topGenes_manyAsPercent_g_rootIdx) == 1)
    # topGenes_manyAsPercent_g_rootGO <- names(V(topGenes_manyAsPercent_g)[topGenes_manyAsPercent_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topGenes_manyAsPercent_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topGenes_manyAsPercent_g_density <- edge_density(topGenes_manyAsPercent_g)
    
    # edge connectivity (adhesion)
    # topGenes_manyAsPercent_g_undir <- as.undirected(topGenes_manyAsPercent_g)
    topGenes_manyAsPercent_g_adhesion <- edge_connectivity(topGenes_manyAsPercent_g)
    
    # vertex connectivity (cohesion)
    topGenes_manyAsPercent_g_cohesion <- vertex_connectivity(topGenes_manyAsPercent_g)
    # diameter
    topGenes_manyAsPercent_g_diameter <- diameter(topGenes_manyAsPercent_g)
    
    # mean geodesic distance
    topGenes_manyAsPercent_g_meanDist <- mean_distance(topGenes_manyAsPercent_g, directed = FALSE) # directed TRUE or FALSE ???
    topGenes_manyAsPercent_g_meanDistDir <- mean_distance(topGenes_manyAsPercent_g, directed = TRUE) # directed TRUE or FALSE ???
    
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topGenes_manyAsPercent_g_meanEccentricity <- mean(eccentricity(topGenes_manyAsPercent_g))

   # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topGenes_manyAsPercent_g_meanBetweenness <- mean(betweenness(topGenes_manyAsPercent_g))
  
      
  } else {
    topGenes_manyAsPercent_signifGOmin <- NA
    nTopGenes_manyAsPercent_signifGOmin <- NA
    topGenes_manyAsPercent_g_nNodes <- NA
    topGenes_manyAsPercent_g_density <- NA
    topGenes_manyAsPercent_g_adhesion <- NA
    topGenes_manyAsPercent_g_cohesion <- NA
    topGenes_manyAsPercent_g_diameter <- NA
    topGenes_manyAsPercent_g_meanDist <- NA
    topGenes_manyAsPercent_g_meanDistDir <- NA
    topGenes_manyAsPercent_g_meanEccentricity <- NA
    topGenes_manyAsPercent_g_meanBetweenness <- NA
  }
    
    
  ############################################################################################
  ############################################################################################ GRAPH - manyAsTopTADs
  ############################################################################################
  
  if(!is.na(nTopGenes_manyAsTopTADs_signifGOid) & nTopGenes_manyAsTopTADs_signifGOid > 0) {
    
    topGenes_manyAsTopTADs_signifGOmin <- minimal_set(go, topGenes_manyAsTopTADs_signifGOid)
    nTopGenes_manyAsTopTADs_signifGOmin <- length(topGenes_manyAsTopTADs_signifGOmin)
    topGenes_manyAsTopTADs_gNEL <- inducedGraph(dag = GO_g, startNodes = topGenes_manyAsTopTADs_signifGOmin[topGenes_manyAsTopTADs_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topGenes_manyAsTopTADs_gNEL_rev <- reverseArch(topGenes_manyAsTopTADs_gNEL) # topGO
    topGenes_manyAsTopTADs_gNEL_rev_noAll <- removeNode("all", topGenes_manyAsTopTADs_gNEL_rev) # not needed ? -> do not want this root
    topGenes_manyAsTopTADs_g <- igraph.from.graphNEL(graphNEL = topGenes_manyAsTopTADs_gNEL_rev_noAll)
    
    stopifnot(is.connected(topGenes_manyAsTopTADs_g))
    stopifnot(is.directed(topGenes_manyAsTopTADs_g))
    stopifnot(is.dag(topGenes_manyAsTopTADs_g))
    
    topGenes_manyAsTopTADs_g_nNodes <- length(V(topGenes_manyAsTopTADs_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topGenes_manyAsTopTADs_g_rootIdx <- which(sapply(sapply(V(topGenes_manyAsTopTADs_g), 
    #                                 function(x) neighbors(topGenes_manyAsTopTADs_g, x, mode="in")), length) == 0)
    # stopifnot(length(topGenes_manyAsTopTADs_g_rootIdx) == 1)
    # topGenes_manyAsTopTADs_g_rootGO <- names(V(topGenes_manyAsTopTADs_g)[topGenes_manyAsTopTADs_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topGenes_manyAsTopTADs_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topGenes_manyAsTopTADs_g_density <- edge_density(topGenes_manyAsTopTADs_g)
    
    # edge connectivity (adhesion)
    # topGenes_manyAsTopTADs_g_undir <- as.undirected(topGenes_manyAsTopTADs_g)
    topGenes_manyAsTopTADs_g_adhesion <- edge_connectivity(topGenes_manyAsTopTADs_g)
    
    # vertex connectivity (cohesion)
    topGenes_manyAsTopTADs_g_cohesion <- vertex_connectivity(topGenes_manyAsTopTADs_g)
    # diameter
    topGenes_manyAsTopTADs_g_diameter <- diameter(topGenes_manyAsTopTADs_g)
    
    # mean geodesic distance
    topGenes_manyAsTopTADs_g_meanDist <- mean_distance(topGenes_manyAsTopTADs_g, directed = FALSE) # directed TRUE or FALSE ???
    topGenes_manyAsTopTADs_g_meanDistDir <- mean_distance(topGenes_manyAsTopTADs_g, directed = TRUE) # directed TRUE or FALSE ???
    
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topGenes_manyAsTopTADs_g_meanEccentricity <- mean(eccentricity(topGenes_manyAsTopTADs_g))

    # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topGenes_manyAsTopTADs_g_meanBetweenness <- mean(betweenness(topGenes_manyAsTopTADs_g))
    
  } else {
    topGenes_manyAsTopTADs_signifGOmin <- NA
    nTopGenes_manyAsTopTADs_signifGOmin <- NA
    topGenes_manyAsTopTADs_g_nNodes <- NA
    topGenes_manyAsTopTADs_g_density <- NA
    topGenes_manyAsTopTADs_g_adhesion <- NA
    topGenes_manyAsTopTADs_g_cohesion <- NA
    topGenes_manyAsTopTADs_g_diameter <- NA
    topGenes_manyAsTopTADs_g_meanDist <- NA
    topGenes_manyAsTopTADs_g_meanDistDir <- NA
    topGenes_manyAsTopTADs_g_meanEccentricity <- NA
    topGenes_manyAsTopTADs_g_meanBetweenness <- NA
  }
    
    
    


    data.frame(
      dataset = curr_ds, 
      
    nTopTADs_genes=nTopTADs_genes,
    nTopTADs_genes_signifGOid=nTopTADs_genes_signifGOid,
    nTopTADs_genes_signifGOterm=nTopTADs_genes_signifGOterm,
    nTopTADs_genes_signifGOmin=nTopTADs_genes_signifGOmin,
    
    topTADs_genes_g_density=topTADs_genes_g_density,
    topTADs_genes_g_adhesion=topTADs_genes_g_adhesion,
    topTADs_genes_g_cohesion=topTADs_genes_g_cohesion,
    topTADs_genes_g_diameter=topTADs_genes_g_diameter,
    topTADs_genes_g_meanDist=topTADs_genes_g_meanDist,
    topTADs_genes_g_meanDistDir=topTADs_genes_g_meanDistDir,
    topTADs_genes_g_meanEccentricity=topTADs_genes_g_meanEccentricity,
    topTADs_genes_g_meanBetweenness=topTADs_genes_g_meanBetweenness,
    
    nTopGenes_manyAsPercent=nTopGenes_manyAsPercent,
    nTopGenes_manyAsPercent_signifGOid=nTopGenes_manyAsPercent_signifGOid,
    nTopGenes_manyAsPercent_signifGOterm=nTopGenes_manyAsPercent_signifGOterm,
    nTopGenes_manyAsPercent_signifGOmin=nTopGenes_manyAsPercent_signifGOmin,
    
    topGenes_manyAsPercent_g_density=topGenes_manyAsPercent_g_density,
    topGenes_manyAsPercent_g_adhesion=topGenes_manyAsPercent_g_adhesion,
    topGenes_manyAsPercent_g_cohesion=topGenes_manyAsPercent_g_cohesion,
    topGenes_manyAsPercent_g_diameter=topGenes_manyAsPercent_g_diameter,
    topGenes_manyAsPercent_g_meanDist=topGenes_manyAsPercent_g_meanDist,
    topGenes_manyAsPercent_g_meanDistDir=topGenes_manyAsPercent_g_meanDistDir,
    topGenes_manyAsPercent_g_meanEccentricity=topGenes_manyAsPercent_g_meanEccentricity,
    topGenes_manyAsPercent_g_meanBetweenness=topGenes_manyAsPercent_g_meanBetweenness,
    
    nTopGenes_manyAsTopTADs=nTopGenes_manyAsTopTADs,
    nTopGenes_manyAsTopTADs_signifGOid=nTopGenes_manyAsTopTADs_signifGOid,
    nTopGenes_manyAsTopTADs_signifGOterm=nTopGenes_manyAsTopTADs_signifGOterm,
    nTopGenes_manyAsTopTADs_signifGOmin=nTopGenes_manyAsTopTADs_signifGOmin,
    
    topGenes_manyAsTopTADs_g_density=topGenes_manyAsTopTADs_g_density,
    topGenes_manyAsTopTADs_g_adhesion=topGenes_manyAsTopTADs_g_adhesion,
    topGenes_manyAsTopTADs_g_cohesion=topGenes_manyAsTopTADs_g_cohesion,
    topGenes_manyAsTopTADs_g_diameter=topGenes_manyAsTopTADs_g_diameter,
    topGenes_manyAsTopTADs_g_meanDist=topGenes_manyAsTopTADs_g_meanDist,
    topGenes_manyAsTopTADs_g_meanDistDir=topGenes_manyAsTopTADs_g_meanDistDir,
    topGenes_manyAsTopTADs_g_meanEccentricity=topGenes_manyAsTopTADs_g_meanEccentricity,
    topGenes_manyAsTopTADs_g_meanBetweenness=topGenes_manyAsTopTADs_g_meanBetweenness,
    
    nIntersectSignifGOid_manyAsPercent = length(intersect(topTADs_genes_signifGOid, topGenes_manyAsPercent_signifGOid)),
    nIntersectSignifGOterm_manyAsPercent = length(intersect(topTADs_genes_signifGOterm, topGenes_manyAsPercent_signifGOterm)),
    nIntersectSignifGOmin_manyAsPercent = length(intersect(topTADs_genes_signifGOmin, topGenes_manyAsPercent_signifGOmin)),
    
    nIntersectSignifGOid_manyAsTopTADs = length(intersect(topTADs_genes_signifGOid, topGenes_manyAsTopTADs_signifGOid)),
    nIntersectSignifGOterm_manyAsTopTADs = length(intersect(topTADs_genes_signifGOterm, topGenes_manyAsTopTADs_signifGOterm)),
    nIntersectSignifGOmin_manyAsTopTADs = length(intersect(topTADs_genes_signifGOmin, topGenes_manyAsTopTADs_signifGOmin)),
    
    stringsAsFactors=FALSE
    
    )
# density
    # The density of a graph is the ratio of the number of edges and the number of possible edges. 
# => if targeted -> low density ???

# connectivity
# The vertex connectivity of a graph or two vertices, this is recently also called group cohesion. 
# The vertex connectivity of a graph is the minimum vertex connectivity of all (ordered) pairs of vertices in the graph. In other words this is the minimum number of vertices needed 
    # to remove to make the graph not strongly connected. (If the graph is not strongly connected then this is zero.) 
# The cohesion of a graph (as defined by White and Harary, see references), is the vertex connectivity of the graph. This is calculated by cohesion.
# 
# These three functions essentially calculate the same measure(s), more precisely vertex_connectivity is the most general, the other two are included only for the ease of using more descriptive function names. 7
# 
    
    
# diameter of a graph is the length of the longest geodesic
    
# The edge connectivity of a pair of vertices (source and target) is the minimum number of edges needed to remove to eliminate all 
# (directed) paths from source to target. edge_connectivity calculates this quantity if both the source and target arguments are
# given (and not NULL).
# 
# The edge connectivity of a graph is the minimum of the edge connectivity of every (ordered) pair of vertices in the graph. 
# edge_connectivity calculates this quantity if neither the source nor the target arguments are given (ie. they are both NULL).
# 
# A set of edge disjoint paths between two vertices is a set of paths between them containing no common edges. 
# The maximum number of edge disjoint paths between two vertices is the same as their edge connectivity.
# 
# The adhesion of a graph is the minimum number of edges needed to remove to obtain a graph which is not strongly connected. 
#    This is the same as the edge connectivity of the graph.
# 
# The three functions documented on this page calculate similar properties, more precisely the most general is edge_connectivity, 
# the others are included only for having more descriptive function names. 

# COHESIVENESS OF BLOCKS IN SOCIAL NETWORKS (White and Harary )
# What we call cohesion is the contribution made by (adding or subtracting) individual members of a group, 
# together with their ties, to holding it together.
# What we call adhesion (edge cohesion) is the contribution   made   to   holding   a   group   together—keeping   membership
# constant—by (adding or subtracting) ties between its members.
                             # )
    #save(enricher_all_results, file = outFile)
    #cat(paste0("... written: ", outFile, "\n"))
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)

  # all_ds_DT <- all_ds_DT[order(all_ds_DT$nIntersectSignifGO, all_ds_DT$nIntersectSignifGOmin),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

# stop("-- ok -- \n")

# [1] "dataset"                        "nTopTADs_genes"             
# [3] "nTopTADs_genes_signifGOid"   "nTopTADs_genes_signifGOterm"
# [5] "nTopTADs_genes_signifGOmin"  "topTADs_genes_g_density"    
# [7] "topTADs_genes_g_adhesion"    "topTADs_genes_g_cohesion"   
# [9] "topTADs_genes_g_diameter"    "topTADs_genes_g_meanDist"   
# [11] "topTADs_genes_g_meanDistDir" "nTopGenes_manyAsPercent"                  
# [13] "nTopGenes_manyAsPercent_signifGOid"        "nTopGenes_manyAsPercent_signifGOterm"     
# [15] "nTopGenes_manyAsPercent_signifGOmin"       "topGenes_manyAsPercent_g_density"         
# [17] "topGenes_manyAsPercent_g_adhesion"         "topGenes_manyAsPercent_g_cohesion"        
# [19] "topGenes_manyAsPercent_g_diameter"         "topGenes_manyAsPercent_g_meanDist"        
# [21] "topGenes_manyAsPercent_g_meanDistDir"      "nIntersectSignifGOid"          
# [23] "nIntersectSignifGOterm"         "nIntersectSignifGOmin"         

stopifnot(!duplicated(all_ds_DT$dataset))  

all_ds_DT$nTopTADs_genes_signifGOidRatio <- all_ds_DT$nTopTADs_genes_signifGOid/all_ds_DT$nTopTADs_genes
all_ds_DT$nTopTADs_genes_signifGOtermRatio <- all_ds_DT$nTopTADs_genes_signifGOterm/all_ds_DT$nTopTADs_genes
all_ds_DT$nTopTADs_genes_signifGOminRatio <- all_ds_DT$nTopTADs_genes_signifGOmin/all_ds_DT$nTopTADs_genes

all_ds_DT$nTopGenes_manyAsPercent_signifGOidRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOid/all_ds_DT$nTopGenes_manyAsPercent
all_ds_DT$nTopGenes_manyAsPercent_signifGOtermRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOterm/all_ds_DT$nTopGenes_manyAsPercent
all_ds_DT$nTopGenes_manyAsPercent_signifGOminRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOmin/all_ds_DT$nTopGenes_manyAsPercent

all_ds_DT$nTopGenes_manyAsTopTADs_signifGOidRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOid/all_ds_DT$nTopGenes_manyAsTopTADs
all_ds_DT$nTopGenes_manyAsTopTADs_signifGOtermRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOterm/all_ds_DT$nTopGenes_manyAsTopTADs
all_ds_DT$nTopGenes_manyAsTopTADs_signifGOminRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOmin/all_ds_DT$nTopGenes_manyAsTopTADs

##########################################################################################

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

myGenSub <- paste0("topPercent = ", topPercent, "; pvalSelectGO = ", pvalSelectGO)

outDT <- all_ds_DT
outDT <- outDT[order(outDT$aucFCC, decreasing=TRUE),]
write.table(outDT, file = file.path(outFold, "all_ds_DT.txt"), sep="\t", quote=F, col.names=T, row.names=F)
#stop("..ok..\n")


##########################################################################################

# scatterplot vs. aucFCC and aucCoexprDist


all_vars <- colnames(all_ds_DT)
all_vars <- all_vars[all_vars != "dataset"]

ref_vars <- c("aucFCC", "aucCoexprDist")

ref_var=ref_vars[1]
curr_var=all_vars[1]

for(ref_var in ref_vars){
  for(curr_var in all_vars) {
    
    if(ref_var == curr_var) next
    
    if(grepl("topTADs_genes", curr_var) | grepl("TopTADs_genes", curr_var)) {
      mysub <- paste0("topTADs_genes - ", myGenSub)
    } else if(grepl("topGenes_manyAsPercent", curr_var) | grepl("TopGenes_manyAsPercent", curr_var)) {
      mysub <- paste0("topGenes_manyAsPercent - ", myGenSub)
    } else if(grepl("topGenes_manyAsTopTADs", curr_var) | grepl("TopGenes_manyAsTopTADs", curr_var)) {
      mysub <- paste0("topGenes_manyAsTopTADs - ", myGenSub)
    }
    
    myylab <- unique(gsub("topTADs_genes_", "", 
                            gsub("topGenes_manyAsPercent_", "", 
                                 gsub("topGenes_manyAsTopTADs_", "", 
                                 gsub("nTopGenes_manyAsPercent_", "",
                                      gsub("nTopGenes_manyAsTopTADs_", "",
                                      gsub("nTopTADs_genes_", "",  curr_var)))))))
    
    mymain <- paste0(myylab, " vs. ", ref_var)
    
    myx <- all_ds_DT[,ref_var]
    myy <- all_ds_DT[,curr_var]
    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x=myx,
         y=myy,
         main = paste0(myylab, " vs. ", ref_var),
         xlab =  paste0(ref_var),
         ylab = myylab,
         pch=16,cex=0.7)
    add_curv_fit(x=myx,
                 y=myy, withR2 = F, lty=2)
    addCorr(x=myx,
            y=myy, legPos = "topleft", bty="n")
    text(x = myx, y = myy, labels = all_ds_DT[,"dataset"], cex=0.7)
    mtext(mysub,side=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

########################################## ########################################## 

all_cmp_var <- c("topGenes_manyAsPercent",  "topGenes_manyAsTopTADs", "topTADs_genes")

all_vars_toplot <- unique(gsub("topTADs_genes_", "", 
                      gsub("topGenes_manyAsPercent_", "", 
                           gsub("topGenes_manyAsTopTADs_", "", 
                           gsub("nTopGenes_manyAsPercent_", "",
                                gsub("nTopGenes_manyAsTopTADs_", "",
                                gsub("nTopTADs_genes_", "",  colnames(all_ds_DT))))))))

all_vars_toplot <- all_vars_toplot[! all_vars_toplot %in% c("dataset", "aucFCC", "aucCoexprDist",
                                                            "nTopGenes_manyAsPercent",  "nTopGenes_manyAsTopTADs", 
                                                            "nTopTADs_genes",
                                                            "nIntersectSignifGOterm_manyAsPercent","nIntersectSignifGOterm_manyAsTopTADs",
                                                            "nIntersectSignifGOid_manyAsPercent", "nIntersectSignifGOid_manyAsTopTADs", 
                                                            "nIntersectSignifGOmin_manyAsPercent","nIntersectSignifGOmin_manyAsTopTADs") ]

all_vars_toplot <- c(all_vars_toplot, "nTop")

var_to_plot=all_vars_toplot[10]
cmp_var=all_cmp_var[10]

for(var_to_plot in all_vars_toplot){
  
  if(var_to_plot == "nTop") {
    cat("A\n")
    # plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, ".+enes_"), colnames(all_ds_DT))]    
    plotDT <- all_ds_DT[, grepl("^nTopTADs_genes$|^nTopGenes_manyAsTopTADs$|^nTopGenes_manyAsPercent$", colnames(all_ds_DT))]
  } else {
    # cat("B\n")
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, "$"), colnames(all_ds_DT))]
  }

  plotDT$dataset <- all_ds_DT$dataset
  plotDT_m <- melt(plotDT, id="dataset")
  
  plotDT_m$dataset <- factor(plotDT_m$dataset, levels = all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=TRUE)])
  
  plotDT_m$variable_leg <- gsub(paste0("_", var_to_plot), "", plotDT_m$variable)
  
  plotTit <- paste0(var_to_plot)
  mySub <- paste0(myGenSub)
    
  # stopifnot(length(unique(plotDT_m$variable)) == 2)
  stopifnot(length(unique(plotDT_m$variable)) == 3)
    
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(plotDT_m$variable)),
                        # values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$variable)),
                        labels = setNames(unique(plotDT_m$variable_leg), unique(plotDT_m$variable)))+
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
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_barplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    plotDT_m$variable_leg <- as.character(plotDT_m$variable_leg)
    
#    my_comparisons <- list(unique(plotDT_m$variable_leg))

x <- combn(unique(plotDT_m$variable_leg),2)
my_comparisons <- lapply(seq_len(ncol(x)), function(i) unlist(x[,i]))

    p_box <- ggviolin(plotDT_m, x = "variable_leg", y = "value", 
                      fill = "variable_leg",
                      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                      title = plotTit,
                      xlab="",
                      ylab = paste0(var_to_plot),# gsub("_", " ", var_to_plot),
                      # sub=mySub,
                      legend.title = "",
                      add = "boxplot", add.params = list(fill = "white"))+
      stat_compare_means(comparisons = my_comparisons, 
                         method="wilcox",
                         # position = "bottomleft",
                         label = "p.format")    + # Add significance levels
      # stat_compare_means()
    
      # stat_compare_means(comparisons = my_comparisons, 
      #                    # position = "bottomleft",
      #                    label = c("p.format"))+ # Add significance levels
      # stat_compare_means(comparisons = my_comparisons,
      #                    # position="bottom",
      #                    label.x.npc = 0,
      #                    label.y.npc=0,
      #                    method="wilcox",
      #                    label = c("p.format"))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot.", plotType))
    ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
}



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

