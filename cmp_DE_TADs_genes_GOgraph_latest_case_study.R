startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GOgraph_latest_case_study.R\n"))

# Rscript cmp_DE_TADs_genes_GOgraph_latest_case_study.R

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
source("GO_graph_utils.R")

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
topPercent <- 0.10
pvalSelect <- 0.05
pvalSelectGO <- 0.05

curr_ds <- "TCGAcrc_msi_mss"

enricher_ontologyType <- "BP"

outFold <- file.path("CMP_TADs_GENES_GOgraph_latest_CASE_STUDY",  enricher_ontologyType, paste0(topPercent, "_", pvalSelect, "_", pvalSelectGO), curr_ds)
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

stopifnot(is.connected(selectTADs_genes_g))
stopifnot(is.directed(selectTADs_genes_g))
stopifnot(is.dag(selectTADs_genes_g))

selectTADs_genes_g_nNodes <- length(V(selectTADs_genes_g))

# density
selectTADs_genes_g_density <- edge_density(selectTADs_genes_g)

# edge connectivity (adhesion)
# selectTADs_genes_g_undir <- as.undirected(selectTADs_genes_g)
selectTADs_genes_g_adhesion <- edge_connectivity(selectTADs_genes_g)

# vertex connectivity (cohesion)
selectTADs_genes_g_cohesion <- vertex_connectivity(selectTADs_genes_g)
# diameter
selectTADs_genes_g_diameter <- diameter(selectTADs_genes_g)

# mean geodesic distance
selectTADs_genes_g_meanDist <- mean_distance(selectTADs_genes_g, directed = FALSE) # directed TRUE or FALSE ???
selectTADs_genes_g_meanDistDir <- mean_distance(selectTADs_genes_g, directed = TRUE) # directed TRUE or FALSE ???



############################################################################################
############################################################################################ selectGenes
############################################################################################
#***** 2) selectGenes

selectGenes_enrich <- enricher(gene = selectGenes, 
                             TERM2GENE=c5_msigdb,
                             pvalueCutoff = enricher_pvalueCutoff, 
                             pAdjustMethod = enricher_pAdjustMethod, 
                             minGSSize = enricher_minGSSize, 
                             maxGSSize = enricher_maxGSSize, 
                             qvalueCutoff =enricher_qvalueCutoff)

selectGenes_resultDT <- selectGenes_enrich@result
selectGenes_resultDT <- selectGenes_resultDT[order(selectGenes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

selectGenes_signifGOterm <- as.character(selectGenes_resultDT$ID[selectGenes_resultDT$p.adjust <= pvalSelectGO])
nSelectGenes_signifGOterm <- length(selectGenes_signifGOterm)

selectGenes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", selectGenes_signifGOterm)))

length(selectGenes_signifGOterm_tmp)
# 5
sum(selectGenes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
# 5
sum(selectGenes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
# 5
sum(selectGenes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
# 5

selectGenes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% selectGenes_signifGOterm_tmp]
nSelectGenes_signifGOid <- length(selectGenes_signifGOid)

###### GO analysis

selectGenes_signifGOmin <- minimal_set(go, selectGenes_signifGOid)
nSelectGenes_signifGOmin <- length(selectGenes_signifGOmin)
selectGenes_gNEL <- inducedGraph(dag = GO_g, startNodes = selectGenes_signifGOmin[selectGenes_signifGOmin %in% nodes(GO_g)])

# !!! need to reverse the graph !!!
selectGenes_gNEL_rev <- reverseArch(selectGenes_gNEL) # topGO
selectGenes_gNEL_rev_noAll <- removeNode("all", selectGenes_gNEL_rev) # not needed ? -> do not want this root
selectGenes_g <- igraph.from.graphNEL(graphNEL = selectGenes_gNEL_rev_noAll)

stopifnot(is.connected(selectGenes_g))
stopifnot(is.directed(selectGenes_g))
stopifnot(is.dag(selectGenes_g))

selectGenes_g_nNodes <- length(V(selectGenes_g))

# density
selectGenes_g_density <- edge_density(selectGenes_g)

# edge connectivity (adhesion)
# selectGenes_g_undir <- as.undirected(selectGenes_g)
selectGenes_g_adhesion <- edge_connectivity(selectGenes_g)

# vertex connectivity (cohesion)
selectGenes_g_cohesion <- vertex_connectivity(selectGenes_g)
# diameter
selectGenes_g_diameter <- diameter(selectGenes_g)

# mean geodesic distance
selectGenes_g_meanDist <- mean_distance(selectGenes_g, directed = FALSE) # directed TRUE or FALSE ???
selectGenes_g_meanDistDir <- mean_distance(selectGenes_g, directed = TRUE) # directed TRUE or FALSE ???
  
# cohesive_blocks(as.undirected(selectGenes_g))
# plot_hierarchy(cohesive_blocks(as.undirected(selectGenes_g)))
# count_subgraph_isomorphisms(pattern=selectGenes_g, target=selectTADs_genes_g)
# count_subgraph_isomorphisms(pattern=selectTADs_genes_g, target=selectGenes_g)
                            


                            
############################################################################################
############################################################################################ topTADs
############################################################################################

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

###### GO analysis

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

############################################################################################
############################################################################################ manyAsPercent
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

###### GO analysis

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

############################################################################################
############################################################################################ manyAsTopTADs
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

###### GO analysis

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

# density
  # The density of a graph is the ratio of the number of edges and the number of possible edges. 
# => if targeted -> low density ???

# connectivity
# The vertex connectivity of a graph or two vertices, this is recently also called group cohesion. 
# The vertex connectivity of a graph is the minimum vertex connectivity of all (ordered) pairs of vertices in the graph. In other words this is the minimum number of vertices needed to remove to make the graph not strongly connected. (If the graph is not strongly connected then this is zero.) 
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

#======================================================================================================================
#====================================================================================================================== draw the graph for the top GO
#======================================================================================================================


#********************************************** selectTADs_genes
#### for visualization, select the top GO
nTop <- 1


#********************************************** selectTADs_genes
# associate the Term with pvals
selectTADs_genes_signifGOterm_pvals <- setNames(  selectTADs_genes_resultDT$p.adjust[order(selectTADs_genes_resultDT$p.adjust, decreasing=FALSE)],
                                                  tolower(gsub("_", " ", gsub("^GO_", "", selectTADs_genes_resultDT$ID[order(selectTADs_genes_resultDT$p.adjust, decreasing=FALSE)]))))
# take only those for which I have the ids
selectTADs_genes_signifGOterm_pvals <- selectTADs_genes_signifGOterm_pvals[names(selectTADs_genes_signifGOterm_pvals) %in% names(GOterm_GOid_cP)]
# associate the IDs with pvals
selectTADs_genes_signifGOid_pvals <- selectTADs_genes_signifGOterm_pvals
names(selectTADs_genes_signifGOid_pvals) <- sapply(names(selectTADs_genes_signifGOterm_pvals), function(x)  GOterm_GOid_cP[names(GOterm_GOid_cP) == x])
stopifnot(!any(length(names(selectTADs_genes_signifGOid_pvals)) == 0))
# take the ones that are in the minimal set IDs
selectTADs_genes_signifGOminID_pvals  <- selectTADs_genes_signifGOid_pvals[names(selectTADs_genes_signifGOid_pvals) %in% selectTADs_genes_signifGOmin]
# take the top
selectTADs_genes_topSignifGOmin <- names(selectTADs_genes_signifGOminID_pvals[1:nTop])
stopifnot(length(selectTADs_genes_topSignifGOmin) == 1)
selectTADs_genes_id_to_plot <- selectTADs_genes_topSignifGOmin[selectTADs_genes_topSignifGOmin %in% nodes(GO_g)]
stopifnot(length(selectTADs_genes_id_to_plot) == 1)

outFile <- file.path(outFold, paste0("selectTADs_genes_GOgraph.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
my_plot_GO_graph(GOgraph = GO_g, GOids = selectTADs_genes_id_to_plot)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#********************************************** selectGenes
# associate the Term with pvals
selectGenes_signifGOterm_pvals <- setNames(  selectGenes_resultDT$p.adjust[order(selectGenes_resultDT$p.adjust, decreasing=FALSE)],
                                             tolower(gsub("_", " ", gsub("^GO_", "", selectGenes_resultDT$ID[order(selectGenes_resultDT$p.adjust, decreasing=FALSE)]))))
# take only those for which I have the ids
selectGenes_signifGOterm_pvals <- selectGenes_signifGOterm_pvals[names(selectGenes_signifGOterm_pvals) %in% names(GOterm_GOid_cP)]
# associate the IDs with pvals
selectGenes_signifGOid_pvals <- selectGenes_signifGOterm_pvals
names(selectGenes_signifGOid_pvals) <- sapply(names(selectGenes_signifGOterm_pvals), function(x)  GOterm_GOid_cP[names(GOterm_GOid_cP) == x])
stopifnot(!any(length(names(selectGenes_signifGOid_pvals)) == 0))
# take the ones that are in the minimal set IDs
selectGenes_signifGOminID_pvals  <- selectGenes_signifGOid_pvals[names(selectGenes_signifGOid_pvals) %in% selectGenes_signifGOmin]
# take the top
selectGenes_topSignifGOmin <- names(selectGenes_signifGOminID_pvals[1:nTop])
stopifnot(length(selectGenes_topSignifGOmin) == 1)
selectGenes_id_to_plot <- selectGenes_topSignifGOmin[selectGenes_topSignifGOmin %in% nodes(GO_g)]
stopifnot(length(selectGenes_id_to_plot) == 1)

outFile <- file.path(outFold, paste0("selectGenes_GOgraph.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
my_plot_GO_graph(GOgraph = GO_g, GOids = selectGenes_id_to_plot)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#********************************************** topTADs
# associate the Term with pvals
topTADs_genes_signifGOterm_pvals <- setNames(  topTADs_genes_resultDT$p.adjust[order(topTADs_genes_resultDT$p.adjust, decreasing=FALSE)],
                                               tolower(gsub("_", " ", gsub("^GO_", "", topTADs_genes_resultDT$ID[order(topTADs_genes_resultDT$p.adjust, decreasing=FALSE)]))))
# take only those for which I have the ids
topTADs_genes_signifGOterm_pvals <- topTADs_genes_signifGOterm_pvals[names(topTADs_genes_signifGOterm_pvals) %in% names(GOterm_GOid_cP)]
# associate the IDs with pvals
topTADs_genes_signifGOid_pvals <- topTADs_genes_signifGOterm_pvals
names(topTADs_genes_signifGOid_pvals) <- sapply(names(topTADs_genes_signifGOterm_pvals), function(x)  GOterm_GOid_cP[names(GOterm_GOid_cP) == x])
stopifnot(!any(length(names(topTADs_genes_signifGOid_pvals)) == 0))
# take the ones that are in the minimal set IDs
topTADs_genes_signifGOminID_pvals  <- topTADs_genes_signifGOid_pvals[names(topTADs_genes_signifGOid_pvals) %in% topTADs_genes_signifGOmin]
# take the top
topTADs_genes_topSignifGOmin <- names(topTADs_genes_signifGOminID_pvals[1:nTop])
stopifnot(length(topTADs_genes_topSignifGOmin) == 1)
topTADs_genes_id_to_plot <- topTADs_genes_topSignifGOmin[topTADs_genes_topSignifGOmin %in% nodes(GO_g)]
stopifnot(length(topTADs_genes_id_to_plot) == 1)

outFile <- file.path(outFold, paste0("topTADs_genes_GOgraph.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
my_plot_GO_graph(GOgraph = GO_g, GOids = topTADs_genes_id_to_plot)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#********************************************** topGenes_manyAsPercent

# associate the Term with pvals
topGenes_manyAsPercent_signifGOterm_pvals <- setNames(  topGenes_manyAsPercent_resultDT$p.adjust[order(topGenes_manyAsPercent_resultDT$p.adjust, decreasing=FALSE)],
                                                        tolower(gsub("_", " ", gsub("^GO_", "", topGenes_manyAsPercent_resultDT$ID[order(topGenes_manyAsPercent_resultDT$p.adjust, decreasing=FALSE)]))))
# take only those for which I have the ids
topGenes_manyAsPercent_signifGOterm_pvals <- topGenes_manyAsPercent_signifGOterm_pvals[names(topGenes_manyAsPercent_signifGOterm_pvals) %in% names(GOterm_GOid_cP)]
# associate the IDs with pvals
topGenes_manyAsPercent_signifGOid_pvals <- topGenes_manyAsPercent_signifGOterm_pvals
names(topGenes_manyAsPercent_signifGOid_pvals) <- sapply(names(topGenes_manyAsPercent_signifGOterm_pvals), function(x)  GOterm_GOid_cP[names(GOterm_GOid_cP) == x])
stopifnot(!any(length(names(topGenes_manyAsPercent_signifGOid_pvals)) == 0))
# take the ones that are in the minimal set IDs
topGenes_manyAsPercent_signifGOminID_pvals  <- topGenes_manyAsPercent_signifGOid_pvals[names(topGenes_manyAsPercent_signifGOid_pvals) %in% topGenes_manyAsPercent_signifGOmin]
# take the top
topGenes_manyAsPercent_topSignifGOmin <- names(topGenes_manyAsPercent_signifGOminID_pvals[1:nTop])
stopifnot(length(topGenes_manyAsPercent_topSignifGOmin) == 1)
topGenes_manyAsPercent_id_to_plot <- topGenes_manyAsPercent_topSignifGOmin[topGenes_manyAsPercent_topSignifGOmin %in% nodes(GO_g)]
stopifnot(length(topGenes_manyAsPercent_id_to_plot) == 1)

outFile <- file.path(outFold, paste0("topGenes_manyAsPercent_GOgraph.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
my_plot_GO_graph(GOgraph = GO_g, GOids = topGenes_manyAsPercent_id_to_plot)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#********************************************** topGenes_manyAsTopTADs

# associate the Term with pvals
topGenes_manyAsTopTADs_signifGOterm_pvals <- setNames(  topGenes_manyAsTopTADs_resultDT$p.adjust[order(topGenes_manyAsTopTADs_resultDT$p.adjust, decreasing=FALSE)],
                                                        tolower(gsub("_", " ", gsub("^GO_", "", topGenes_manyAsTopTADs_resultDT$ID[order(topGenes_manyAsTopTADs_resultDT$p.adjust, decreasing=FALSE)]))))
# take only those for which I have the ids
topGenes_manyAsTopTADs_signifGOterm_pvals <- topGenes_manyAsTopTADs_signifGOterm_pvals[names(topGenes_manyAsTopTADs_signifGOterm_pvals) %in% names(GOterm_GOid_cP)]
# associate the IDs with pvals
topGenes_manyAsTopTADs_signifGOid_pvals <- topGenes_manyAsTopTADs_signifGOterm_pvals
names(topGenes_manyAsTopTADs_signifGOid_pvals) <- sapply(names(topGenes_manyAsTopTADs_signifGOterm_pvals), function(x)  GOterm_GOid_cP[names(GOterm_GOid_cP) == x])
stopifnot(!any(length(names(topGenes_manyAsTopTADs_signifGOid_pvals)) == 0))
# take the ones that are in the minimal set IDs
topGenes_manyAsTopTADs_signifGOminID_pvals  <- topGenes_manyAsTopTADs_signifGOid_pvals[names(topGenes_manyAsTopTADs_signifGOid_pvals) %in% topGenes_manyAsTopTADs_signifGOmin]
# take the top
topGenes_manyAsTopTADs_topSignifGOmin <- names(topGenes_manyAsTopTADs_signifGOminID_pvals[1:nTop])
stopifnot(length(topGenes_manyAsTopTADs_topSignifGOmin) == 1)
topGenes_manyAsTopTADs_id_to_plot <- topGenes_manyAsTopTADs_topSignifGOmin[topGenes_manyAsTopTADs_topSignifGOmin %in% nodes(GO_g)]
stopifnot(length(topGenes_manyAsTopTADs_id_to_plot) == 1)

outFile <- file.path(outFold, paste0("topGenes_manyAsTopTADs_GOgraph.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
my_plot_GO_graph(GOgraph = GO_g, GOids = topGenes_manyAsTopTADs_id_to_plot)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

