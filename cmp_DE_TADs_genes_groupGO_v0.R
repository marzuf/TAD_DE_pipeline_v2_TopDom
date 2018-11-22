startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_groupGO.R\n"))

# Rscript cmp_DE_TADs_genes_groupGO_pvalSelect.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- TRUE

library(reshape2)
library(ggplot2)


plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

source("analysis_utils.R")


topPercent <- 0.05
# cmp the top topPercent TADs and topPercent genes


outFold <- "CMP_TADs_GENES_groupGO"
if(!SSHFS) system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_groupgo_logFile.txt")
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

semMeasure <- "Rel"
combMeth <- "avg"
sim_orgDB <- "org.Hs.eg.db"
sim_keytype <- "ENTREZID"
groupGO_ont <- "BP"
groupGO_level <- 3

# needed for GOSemSim
hsGO <- godata('org.Hs.eg.db', ont=groupGO_ont)

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... semMeasure:\t", semMeasure, "\n")
printAndLog(txt, logFile)
txt <- paste0("... combMeth:\t", combMeth, "\n")
printAndLog(txt, logFile)
txt <- paste0("... sim_orgDB:\t", sim_orgDB, "\n")
printAndLog(txt, logFile)
txt <- paste0("... sim_keytype:\t", sim_keytype, "\n")
printAndLog(txt, logFile)
txt <- paste0("... groupGO_ont:\t", groupGO_ont, "\n")
printAndLog(txt, logFile)
txt <- paste0("... groupGO_level:\t", groupGO_level, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:10]
all_ds <- all_ds[1:2]
curr_ds = all_ds[1]
all_ds=c("GSE66306_before_after", "GSE65540_before_after")
curr_ds="TCGAcrc_msi_mss"
# curr_ds = all_ds[1]
if(buildTable){
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
    



    ####### TOP TADs GENES
      cat("... A) start topTADs_genes semantic similarity \n")
      # return similarity value
      topTADs_genes_geneSim <- mgeneSim(genes=topTADs_genes,
                                           semData=hsGO,
                                           measure=semMeasure,
                                           combine = combMeth,
                                           verbose=FALSE)
      
      cat("... B) start groupGO topTADs_genes\n")
      topTADs_genes_groupGO <- groupGO(gene     = topTADs_genes,
                                          OrgDb    = sim_orgDB,
                                          keytype = sim_keytype,
                                          ont      = groupGO_ont,
                                          level    = groupGO_level,
                                          readable = FALSE)
      topTADs_genes_groupGO_nGO <- length(unique(topTADs_genes_groupGO@result$Description))
      topTADs_genes_groupGO_GO <- as.character(topTADs_genes_groupGO@result$ID)
      stopifnot(!duplicated(as.character(topTADs_genes_groupGO@result$ID)))
      
      cat("... C) start topTADs_genes_groupGO_GO semantic similarity \n")
      topTADs_genes_groupGO_GO_sim_mat <- mgoSim(GO1=topTADs_genes_groupGO_GO, 
                                                    GO2=topTADs_genes_groupGO_GO, 
                                                    semData=hsGO, 
                                                    measure=semMeasure, 
                                                    combine=NULL)
      stopifnot(colnames(topTADs_genes_groupGO_GO_sim_mat) == rownames(topTADs_genes_groupGO_GO_sim_mat))
      topTADs_genes_groupGO_GO_sim <- mean(topTADs_genes_groupGO_GO_sim_mat[lower.tri(topTADs_genes_groupGO_GO_sim_mat, diag=F)])



    ####### TOP GENES MANY AS TOP TADs
      cat("... D) start topGenes_manyAsTopTADs semantic similarity \n")
      # return similarity value
      topGenes_manyAsTopTADs_geneSim <- mgeneSim(genes=topGenes_manyAsTopTADs,
                                           semData=hsGO,
                                           measure=semMeasure,
                                           combine = combMeth,
                                           verbose=FALSE)
      
      cat("... E) start groupGO topGenes_manyAsTopTADs\n")
      topGenes_manyAsTopTADs_groupGO <- groupGO(gene     = topGenes_manyAsTopTADs,
                                          OrgDb    = sim_orgDB,
                                          keytype = sim_keytype,
                                          ont      = groupGO_ont,
                                          level    = groupGO_level,
                                          readable = FALSE)
      topGenes_manyAsTopTADs_groupGO_nGO <- length(unique(topGenes_manyAsTopTADs_groupGO@result$Description))
      topGenes_manyAsTopTADs_groupGO_GO <- as.character(topGenes_manyAsTopTADs_groupGO@result$ID)
      stopifnot(!duplicated(as.character(topGenes_manyAsTopTADs_groupGO@result$ID)))
      
      cat("... F) start topGenes_manyAsTopTADs_groupGO_GO semantic similarity \n")
      topGenes_manyAsTopTADs_groupGO_GO_sim_mat <- mgoSim(GO1=topGenes_manyAsTopTADs_groupGO_GO, 
                                                    GO2=topGenes_manyAsTopTADs_groupGO_GO, 
                                                    semData=hsGO, 
                                                    measure=semMeasure, 
                                                    combine=NULL)
      stopifnot(colnames(topGenes_manyAsTopTADs_groupGO_GO_sim_mat) == rownames(topGenes_manyAsTopTADs_groupGO_GO_sim_mat))
      topGenes_manyAsTopTADs_groupGO_GO_sim <- mean(topGenes_manyAsTopTADs_groupGO_GO_sim_mat[lower.tri(topGenes_manyAsTopTADs_groupGO_GO_sim_mat, diag=F)])


    ####### TOP GENES MANY AS PERCENT
      cat("... D) start topGenes_manyAsPercent semantic similarity \n")
      # return similarity value
      topGenes_manyAsPercent_geneSim <- mgeneSim(genes=topGenes_manyAsPercent,
                                           semData=hsGO,
                                           measure=semMeasure,
                                           combine = combMeth,
                                           verbose=FALSE)
      
      cat("... E) start groupGO topGenes_manyAsPercent\n")
      topGenes_manyAsPercent_groupGO <- groupGO(gene     = topGenes_manyAsPercent,
                                          OrgDb    = sim_orgDB,
                                          keytype = sim_keytype,
                                          ont      = groupGO_ont,
                                          level    = groupGO_level,
                                          readable = FALSE)
      topGenes_manyAsPercent_groupGO_nGO <- length(unique(topGenes_manyAsPercent_groupGO@result$Description))
      topGenes_manyAsPercent_groupGO_GO <- as.character(topGenes_manyAsPercent_groupGO@result$ID)
      stopifnot(!duplicated(as.character(topGenes_manyAsPercent_groupGO@result$ID)))
      
      cat("... F) start topGenes_manyAsPercent_groupGO_GO semantic similarity \n")
      topGenes_manyAsPercent_groupGO_GO_sim_mat <- mgoSim(GO1=topGenes_manyAsPercent_groupGO_GO, 
                                                    GO2=topGenes_manyAsPercent_groupGO_GO, 
                                                    semData=hsGO, 
                                                    measure=semMeasure, 
                                                    combine=NULL)
      stopifnot(colnames(topGenes_manyAsPercent_groupGO_GO_sim_mat) == rownames(topGenes_manyAsPercent_groupGO_GO_sim_mat))
      topGenes_manyAsPercent_groupGO_GO_sim <- mean(topGenes_manyAsPercent_groupGO_GO_sim_mat[lower.tri(topGenes_manyAsPercent_groupGO_GO_sim_mat, diag=F)])


    
    if(ntopGenes_manyAsTopTADs > 0 & ntopTADs_genes > 0){
      
      cat("... G) compare topGenes_manyAsTopTADs - topTADs_genes\n")
      
      topTADs_genes_topGenes_manyAsTopTADs_sim <- clusterSim(cluster1 = topGenes_manyAsTopTADs,
                                                     cluster2 = topTADs_genes, 
                                                     semData=hsGO, 
                                                     measure=semMeasure, 
                                                     combine=combMeth)



      cat("... H) compare topGenes_manyAsPercent - topTADs_genes\n")
      
      topTADs_genes_topGenes_manyAsPercent_sim <- clusterSim(cluster1 = topGenes_manyAsPercent,
                                                     cluster2 = topTADs_genes, 
                                                     semData=hsGO, 
                                                     measure=semMeasure, 
                                                     combine=combMeth)




    
    
    data.frame(
      dataset = curr_ds,
      
      ntopTADs = ntopTADs,
      ntopTADs_genes = ntopTADs_genes,
      ntopGenes_manyAsTopTADs = ntopGenes_manyAsTopTADs,
      ntopGenes_manyAsPercent = ntopGenes_manyAsPercent,
      
      topTADs_genes_groupGO_nGO = topTADs_genes_groupGO_nGO,
      topTADs_genes_groupGO_GO_sim = topTADs_genes_groupGO_GO_sim,
      topTADs_genes_geneSim = topTADs_genes_geneSim,
      
      topGenes_manyAsTopTADs_groupGO_nGO = topGenes_manyAsTopTADs_groupGO_nGO,
      topGenes_manyAsTopTADs_groupGO_GO_sim = topGenes_manyAsTopTADs_groupGO_GO_sim,
      topGenes_manyAsTopTADs_geneSim = topGenes_manyAsTopTADs_geneSim,

      topGenes_manyAsPercent_groupGO_nGO = topGenes_manyAsPercent_groupGO_nGO,
      topGenes_manyAsPercent_groupGO_GO_sim = topGenes_manyAsPercent_groupGO_GO_sim,
      topGenes_manyAsPercent_geneSim = topGenes_manyAsPercent_geneSim,
      
      topTADs_genes_topGenes_manyAsTopTADs_sim=topTADs_genes_topGenes_manyAsTopTADs_sim,
      topTADs_genes_topGenes_manyAsPercent_sim=topTADs_genes_topGenes_manyAsPercent_sim,
      
      stringsAsFactors = FALSE
    )
    
    
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  
  orderDataset <- all_ds_DT$dataset[order(all_ds_DT$topTADs_genes_geneSim,
                                          all_ds_DT$topTADs_genes_groupGO_GO_sim,
                                          all_ds_DT$topGenes_manyAsTopTADs_geneSim,
                                          all_ds_DT$topGenes_manyAsTopTADs_groupGO_GO_sim, 
                                          all_ds_DT$topGenes_manyAsPercent_geneSim,
                                          all_ds_DT$topGenes_manyAsPercent_groupGO_GO_sim, decreasing=T),]

  all_ds_DT <- all_ds_DT[match(names(orderDataset), all_ds_DT$dataset),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

stop("...ok...\n")

all_ds_DT$nGO_topTADs_genes_topGenes_manyAsTopTADs <- (all_ds_DT$nbrGO_topTADs_genes + all_ds_DT$nbrGO_topGenes_manyAsTopTADs)
all_ds_DT$ratio_topGenes_manyAsTopTADs <- all_ds_DT$intersectNbrGO_topTADs_genes_topGenes_manyAsTopTADs/all_ds_DT$nGO_topTADs_genes_topGenes_manyAsTopTADs

stopifnot(na.omit(all_ds_DT$ratio_topGenes_manyAsTopTADs) >= 0 & na.omit(all_ds_DT$ratio_topGenes_manyAsTopTADs) <= 1)


########################################## COMPARE % TOP GENES AND SAME # GENES AS TOP TADS
myx <- all_ds_DT$ntopTADs_genes
myy <- all_ds_DT$ntopGenes_manyAsTopTADs
outFile <- file.path(outFold, paste0("cmp_totNbr_diffNbrGenes.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x=myx,
     y=myy,
     main = "# genes compared",
     xlab="genes from DE TADs",
     ylab="DE genes",
     pch=16,cex=0.7)
add_curv_fit(x=myx,
             y=myy, withR2 = F, lty=2)
addCorr(x=myx,
        y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


########################################## BARPLOT RATIO INTERSECT BOTH KIND OF # OF GENES

plotDT <- all_ds_DT[,c("dataset", "ratio_topGenes_manyAsTopTADs")]
plotDT_m <- melt(plotDT, id="dataset")

datasetOrder <- plotDT$dataset[order(plotDT[,2],  decreasing=T)]
plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
stopifnot(!is.na(plotDT_m$dataset))

plotTit <- "Ratio intersect GOs"

p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
  ggtitle(plotTit, subtitle = mySub)+
  geom_bar(stat="identity", position = "dodge")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(""),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values = c(ratio_topGenes_manyAsTopTADs = "dodgerblue4"),
                    labels = c(ratio_topGenes_manyAsTopTADs = "DE genes"))+
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

outFile <- file.path(outFold, paste0("GO_at_intersect_", "bothRatios.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

########################################## BARPLOT # VS. # INTERSECT FOR EACH KIND OF # OF GENES

toplot <- ""

for(curr_plot in toplot) {
    plotDT <- all_ds_DT[,c("dataset", "nGO_topTADs_genes_topGenes_manyAsTopTADs", "intersectNbrGO_topTADs_genes_topGenes_manyAsTopTADs")]
    plotTit <- paste0("# GO and intersect - select DE genes")
    
  plotDT_m <- melt(plotDT, id="dataset")
  
  datasetOrder <- plotDT$dataset[order(plotDT[,3], plotDT[,2], decreasing=T)]
  plotDT_m$dataset <- factor(as.character(plotDT_m$dataset), levels = datasetOrder)
  stopifnot(!is.na(plotDT_m$dataset))
  
  p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
    ggtitle(plotTit, subtitle = mySub)+
    geom_bar(stat="identity", position = "dodge")+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 10))+
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
  
  outFile <- file.path(outFold, paste0("GO_at_intersect_", curr_plot, "TopGenes.", plotType))
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
