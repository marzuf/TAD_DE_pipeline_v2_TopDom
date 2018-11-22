startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO_pvalSelect_latest.R\n"))

# Rscript cmp_DE_TADs_genes_GO_pvalSelect_latest.R

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
pvalSelect <- 0.05
pvalSelectGO <- 0.05

outFold <- "CMP_TADs_GENES_GO_pvalSelect_latest"
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


mySub <- paste0("(pval thresh: signif. genes = ", pvalSelect , "; signif. GO = ",pvalSelectGO, ")") 

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
txt <- paste0("... pval thresh select genes:\t", pvalSelect, "\n")
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
      nSelectTADs_genes <- length(selectTADs_genes)
      
    } else {
      nSelectTADs_genes <- 0
    }
  
  if(nSelectGenes > 0){
    stopifnot(selectGenes %in% gene2tad_DT$entrezID)
    stopifnot(selectGenes %in% pipelineGenes)
  }
    
    # run GO analysis for these 3 sets of entrez genes: selectTADs_genes, selectGenes
    #***** 1) selectTADs_genes
    if(nSelectTADs_genes > 0) {
      selectTADs_genes_enrich <- enricher(gene = selectTADs_genes, 
                                          TERM2GENE=c5_msigdb,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
      
      selectTADs_genes_resultDT <- selectTADs_genes_enrich@result
      selectTADs_genes_resultDT <- selectTADs_genes_resultDT[order(selectTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
      
      # selectTADs_genes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(selectTADs_genes_resultDT)))
      # if(selectTADs_genes_nTop > 0) {
      #   selectTADs_genes_topGo  <- as.character(selectTADs_genes_resultDT$ID[1:selectTADs_genes_nTop])
      # } else {
      #   selectTADs_genes_topGo  <- character(0)  
      # }
      # stopifnot(length(selectTADs_genes_topGo) == selectTADs_genes_nTop)
      
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
      
      # if(nSelectTADs_genes_signifGOterm > 0)
      #   stopifnot(nSelectTADs_genes_signifGOid > 0) # not always TRUE (task 3)
      
    } else {
      selectTADs_genes_resultDT <- data.frame(p.adjust = NA)
      # selectTADs_genes_topGo <- NA
      # selectTADs_genes_nTop <- NA
      selectTADs_genes_signifGOterm <- NA
      selectTADs_genes_signifGOid <- NA
      nSelectTADs_genes_signifGOid <- NA
      nSelectTADs_genes_signifGOterm <- NA
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
      
      # selectGenes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(selectGenes_resultDT)))
      # if(selectGenes_nTop > 0) {
      #   selectGenes_topGo  <- as.character(selectGenes_resultDT$ID[1:selectGenes_nTop])
      # } else {
      #   selectGenes_topGo  <- character(0)  
      # }
      # stopifnot(length(selectGenes_topGo) == selectGenes_nTop)
      
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
      
      # if(nSelectGenes_signifGOterm > 0)
      #   stopifnot(nSelectGenes_signifGOid > 0) # not true task 22
      
    } else {
      
      selectGenes_resultDT <- data.frame(p.adjust = NA)
      # selectGenes_topGo <- NA
      # selectGenes_nTop <- NA
      selectGenes_signifGOid <- NA
      selectGenes_signifGOterm <- NA
      nSelectGenes_signifGOid <- NA
      nSelectGenes_signifGOterm <- NA
    }


###### GO analysis

if(!is.na(nSelectTADs_genes_signifGOid) & nSelectTADs_genes_signifGOid > 0) {
  
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
  
  # Now we can find the root nodes. (either no neighbors or no incident edges)
  # selectTADs_genes_g_rootIdx <- which(sapply(sapply(V(selectTADs_genes_g), 
  #                                 function(x) neighbors(selectTADs_genes_g, x, mode="in")), length) == 0)
  # stopifnot(length(selectTADs_genes_g_rootIdx) == 1)
  # selectTADs_genes_g_rootGO <- names(V(selectTADs_genes_g)[selectTADs_genes_g_rootIdx])
  # txt <- paste0("...... found root GO: ", selectTADs_genes_g_rootGO, "\n")
  # printAndLog(txt, log_file)
  
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

  # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
  selectTADs_genes_g_meanEccentricity <- mean(eccentricity(selectTADs_genes_g))

  # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
  selectTADs_genes_g_meanBetweenness <- mean(betweenness(selectTADs_genes_g))
  
} else {
  selectTADs_genes_signifGOmin <- NA
  nSelectTADs_genes_signifGOmin <- NA
  selectTADs_genes_g_nNodes <- NA
  selectTADs_genes_g_density <- NA
  selectTADs_genes_g_adhesion <- NA
  selectTADs_genes_g_cohesion <- NA
  selectTADs_genes_g_diameter <- NA
  selectTADs_genes_g_meanDist <- NA
  selectTADs_genes_g_meanDistDir <- NA
  selectTADs_genes_g_meanEccentricity <- NA
  selectTADs_genes_g_meanBetweenness <- NA
}
    
    cat("nSelectGenes_signifGOid = ", nSelectGenes_signifGOid, "\n")

if(!is.na(nSelectGenes_signifGOid) & nSelectGenes_signifGOid > 0) {
  
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
  
  # Now we can find the root nodes. (either no neighbors or no incident edges)
  # selectGenes_g_rootIdx <- which(sapply(sapply(V(selectGenes_g), 
  #                                 function(x) neighbors(selectGenes_g, x, mode="in")), length) == 0)
  # stopifnot(length(selectGenes_g_rootIdx) == 1)
  # selectGenes_g_rootGO <- names(V(selectGenes_g)[selectGenes_g_rootIdx])
  # txt <- paste0("...... found root GO: ", selectGenes_g_rootGO, "\n")
  # printAndLog(txt, log_file)
  
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

  # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
  selectGenes_g_meanEccentricity <- mean(eccentricity(selectGenes_g))

  # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
  selectGenes_g_meanBetweenness <- mean(betweenness(selectGenes_g))

    
} else {
  selectGenes_signifGOmin <- NA
  nSelectGenes_signifGOmin <- NA
  selectGenes_g_nNodes <- NA
  selectGenes_g_density <- NA
  selectGenes_g_adhesion <- NA
  selectGenes_g_cohesion <- NA
  selectGenes_g_diameter <- NA
  selectGenes_g_meanDist <- NA
  selectGenes_g_meanDistDir <- NA
  selectGenes_g_meanEccentricity <- NA
  selectGenes_g_meanBetweenness <- NA
}
    
    # if(!is.na(selectTADs_genes_signifGOid) & !is.na(selectTADs_genes_signifGOid))
    # 
    # nIntersectSignifGOid <-  length(intersect(selectTADs_genes_signifGOid, selectGenes_signifGOid)),
    # if(length() == 0) <- 
    # nIntersectSignifGOterm <- length(intersect(selectTADs_genes_signifGOterm, selectGenes_signifGOterm)),
    # if(length() == 0) <- 
    # nIntersectSignifGOmin <- length(intersect(selectTADs_genes_signifGOmin, selectGenes_signifGOmin)),
    # if(length() == 0) <- 
    

    data.frame(
      dataset = curr_ds, 
      
    nSelectTADs_genes=nSelectTADs_genes,
    nSelectTADs_genes_signifGOid=nSelectTADs_genes_signifGOid,
    nSelectTADs_genes_signifGOterm=nSelectTADs_genes_signifGOterm,
    nSelectTADs_genes_signifGOmin=nSelectTADs_genes_signifGOmin,
    
    selectTADs_genes_g_density=selectTADs_genes_g_density,
    selectTADs_genes_g_adhesion=selectTADs_genes_g_adhesion,
    selectTADs_genes_g_cohesion=selectTADs_genes_g_cohesion,
    selectTADs_genes_g_diameter=selectTADs_genes_g_diameter,
    selectTADs_genes_g_meanDist=selectTADs_genes_g_meanDist,
    selectTADs_genes_g_meanDistDir=selectTADs_genes_g_meanDistDir,
    selectTADs_genes_g_meanEccentricity = selectTADs_genes_g_meanEccentricity,
    selectTADs_genes_g_meanBetweenness = selectTADs_genes_g_meanBetweenness,

    nSelectGenes=nSelectGenes,
    nSelectGenes_signifGOid=nSelectGenes_signifGOid,
    nSelectGenes_signifGOterm=nSelectGenes_signifGOterm,
    nSelectGenes_signifGOmin=nSelectGenes_signifGOmin,
    
    selectGenes_g_density=selectGenes_g_density,
    selectGenes_g_adhesion=selectGenes_g_adhesion,
    selectGenes_g_cohesion=selectGenes_g_cohesion,
    selectGenes_g_diameter=selectGenes_g_diameter,
    selectGenes_g_meanDist=selectGenes_g_meanDist,
    selectGenes_g_meanDistDir=selectGenes_g_meanDistDir,
    selectGenes_g_meanEccentricity = selectGenes_g_meanEccentricity,
    selectGenes_g_meanBetweenness = selectGenes_g_meanBetweenness,
    
    nIntersectSignifGOid = length(intersect(selectTADs_genes_signifGOid, selectGenes_signifGOid)),
    nIntersectSignifGOterm = length(intersect(selectTADs_genes_signifGOterm, selectGenes_signifGOterm)),
    nIntersectSignifGOmin = length(intersect(selectTADs_genes_signifGOmin, selectGenes_signifGOmin)),
    
    stringsAsFactors=FALSE
    
    )
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

# [1] "dataset"                        "nSelectTADs_genes"             
# [3] "nSelectTADs_genes_signifGOid"   "nSelectTADs_genes_signifGOterm"
# [5] "nSelectTADs_genes_signifGOmin"  "selectTADs_genes_g_density"    
# [7] "selectTADs_genes_g_adhesion"    "selectTADs_genes_g_cohesion"   
# [9] "selectTADs_genes_g_diameter"    "selectTADs_genes_g_meanDist"   
# [11] "selectTADs_genes_g_meanDistDir" "nSelectGenes"                  
# [13] "nSelectGenes_signifGOid"        "nSelectGenes_signifGOterm"     
# [15] "nSelectGenes_signifGOmin"       "selectGenes_g_density"         
# [17] "selectGenes_g_adhesion"         "selectGenes_g_cohesion"        
# [19] "selectGenes_g_diameter"         "selectGenes_g_meanDist"        
# [21] "selectGenes_g_meanDistDir"      "nIntersectSignifGOid"          
# [23] "nIntersectSignifGOterm"         "nIntersectSignifGOmin"         

stopifnot(!duplicated(all_ds_DT$dataset))

all_ds_DT$nSelectTADs_genes_signifGOidRatio <- all_ds_DT$nSelectTADs_genes_signifGOid/all_ds_DT$nSelectTADs_genes
all_ds_DT$nSelectTADs_genes_signifGOtermRatio <- all_ds_DT$nSelectTADs_genes_signifGOterm/all_ds_DT$nSelectTADs_genes
all_ds_DT$nSelectTADs_genes_signifGOminRatio <- all_ds_DT$nSelectTADs_genes_signifGOmin/all_ds_DT$nSelectTADs_genes

all_ds_DT$nSelectGenes_signifGOidRatio <- all_ds_DT$nSelectGenes_signifGOid/all_ds_DT$nSelectGenes
all_ds_DT$nSelectGenes_signifGOtermRatio <- all_ds_DT$nSelectGenes_signifGOterm/all_ds_DT$nSelectGenes
all_ds_DT$nSelectGenes_signifGOminRatio <- all_ds_DT$nSelectGenes_signifGOmin/all_ds_DT$nSelectGenes


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

myGenSub <- paste0("pvalSelect = ", pvalSelect, "; pvalSelectGO = ", pvalSelectGO)


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
    
    if(grepl("selectTADs_genes", curr_var) | grepl("SelectTADs_genes", curr_var)) {
      mysub <- paste0("selectTADs_genes - ", myGenSub)
    }else if(grepl("selectGenes", curr_var) | grepl("SelectGenes", curr_var)) {
      mysub <- paste0("selectGenes - ", myGenSub)
    }
    
    myylab <- unique(gsub("selectTADs_genes_", "", 
                            gsub("selectGenes_", "", 
                                 gsub("nSelectGenes_", "", 
                                      gsub("nSelectTADs_genes_", "",  curr_var)))))
    
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

all_cmp_var <- c("selectGenes", "selectTADs_genes")

all_vars_toplot <- unique(gsub("selectTADs_genes_", "", 
                      gsub("selectGenes_", "", 
                           gsub("nSelectGenes_", "", 
                                gsub("nSelectTADs_genes_", "",  colnames(all_ds_DT))))))

all_vars_toplot <- all_vars_toplot[! all_vars_toplot %in% c("dataset", "aucFCC", "aucCoexprDist",
                                                            "nSelectGenes", "nSelectTADs_genes",
                                                            "nIntersectSignifGOterm","nIntersectSignifGOid", "nIntersectSignifGOmin" ) ]

all_vars_toplot <- c(all_vars_toplot, "nSelect")

var_to_plot=all_vars_toplot[1]
cmp_var=all_cmp_var[1]

for(var_to_plot in all_vars_toplot){
  
  if(var_to_plot == "nSelect") {
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, ".+enes$"), colnames(all_ds_DT))]    
  } else {
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, "$"), colnames(all_ds_DT))]
  }

  plotDT$dataset <- all_ds_DT$dataset
  plotDT_m <- melt(plotDT, id="dataset")
  
  plotDT_m$dataset <- factor(plotDT_m$dataset, levels = all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=TRUE)])
  
  plotDT_m$variable_leg <- gsub(paste0("_", var_to_plot), "", plotDT_m$variable)
  
  plotTit <- paste0(var_to_plot)
  mySub <- paste0(myGenSub)
    
  stopifnot(length(unique(plotDT_m$variable)) == 2)
    
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$variable)),
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
    
    # my_comparisons <- list(c(ref_var, cmp_var))
    # 
    my_comparisons <- list(unique(plotDT_m$variable_leg))
    # my_comparisons <- list(unique(plotDT_m$variable_leg))
    
    # my_comparisons <- list(c("nSelectTADs_genes", "nSelectGenes"))
    
    # cat(my_comparisons)
    
    p_box <- ggviolin(plotDT_m, x = "variable_leg", y = "value", 
                      fill = "variable_leg",
                      palette = c("#00AFBB", "#FC4E07"),
                      title = plotTit,
                      xlab="",
                      ylab = paste0(var_to_plot),# gsub("_", " ", var_to_plot),
                      # sub=mySub,
                      legend.title = "",
                      add = "boxplot", add.params = list(fill = "white"))+
      # stat_compare_means(comparisons = my_comparisons, 
      #                    # position = "bottomleft",
      #                    label = "p.signif")+ # Add significance levels
      # stat_compare_means()+
      stat_compare_means(comparisons = my_comparisons, 
                         method="wilcox",
                         # position = "bottomleft",
                         label = "p.format")+ # Add significance levels
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

