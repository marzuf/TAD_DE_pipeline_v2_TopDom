startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GOgraph_pvalSelect.R\n"))

# Rscript cmp_DE_TADs_genes_GOgraph_pvalSelect.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology, go, minimal_set
suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(topGO, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  #inducedGraph
suppressPackageStartupMessages(library(AnnotationDbi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # makeGOGraph
suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # makeGOGraph
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # makeGOGraph
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- TRUE

data(go) # ontologyIndex


plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
source("analysis_utils.R")
source("GO_graph_utils.R")

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05

ontologyType <- "MF"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  if(args[1] %in% c("BP", "MF")){
    ontologyType <- args[1]
  }
}

GO_OFFSPRING <- as.list(eval(parse(text = paste0("GO", ontologyType, "OFFSPRING"))))
GO_ANCESTOR <- as.list(eval(parse(text = paste0("GO", ontologyType, "ANCESTOR"))))
GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)


outFold <- file.path("CMP_TADs_GENES_GOgraph_pvalSelect", ontologyType)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_gograph_logFile.txt")
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
txt <- paste0("... pvalSelect:\t", pvalSelect, "\n")
printAndLog(txt, logFile)
txt <- paste0("... ontologyType:\t", ontologyType, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:2]

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
    
    
    ##################
    ###### offAncRatio and tree score -- selectTADs_genes
    ##################
    
    # for each of the selectTADs_genes
    # -> GO
    # --> GOmin
    # ---> for each GO:
    # ----> # all children/(#all children and all parents)    // as a measure of specificity
    # ----> tree score // as measure of targeted
    
    if(nSelectTADs_genes > 0) {
      txt <- paste0("... start selectTADs_genes\n")
      printAndLog(txt, logFile)
      selectTADs_genes_GO_mapping <- bitr(selectTADs_genes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db") # clusterProfiler
      selectTADs_genes_GO_mapping <- selectTADs_genes_GO_mapping[selectTADs_genes_GO_mapping$ONTOLOGY == ontologyType,]
      selectTADs_genes_GO <- unique(selectTADs_genes_GO_mapping$GO)
      selectTADs_genes_GOmin <- minimal_set(go, selectTADs_genes_GO) # ontologyIndex
      
      
      selectTADs_genes_offAncRatio <- get_offAncRatio(GO_list = selectTADs_genes_GOmin, go_data_off=GO_OFFSPRING, go_data_anc=GO_ANCESTOR)
      selectTADs_genes_meanOffAncRatio <- mean(selectTADs_genes_offAncRatio, na.rm=T)
      selectTADs_genes_treeScore <- get_tree_score(GO_list = selectTADs_genes_GOmin, go_graph = GO_g, log_file=logFile)
      
    } else {
      selectTADs_genes_meanOffAncRatio <- NA
      selectTADs_genes_treeScore <- NA
    }

    if(length(selectGenes) > 0) {
      txt <- paste0("... start selectGenes\n")
      printAndLog(txt, logFile)
      selectGenes_GO_mapping <- bitr(selectGenes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db") # clusterProfiler
      selectGenes_GO_mapping <- selectGenes_GO_mapping[selectGenes_GO_mapping$ONTOLOGY == ontologyType,]
      selectGenes_GO <- unique(selectGenes_GO_mapping$GO)
      selectGenes_GOmin <- minimal_set(go, selectGenes_GO) # ontologyIndex
      
      selectGenes_offAncRatio <- get_offAncRatio(GO_list = selectGenes_GOmin, go_data_off=GO_OFFSPRING, go_data_anc=GO_ANCESTOR)
      selectGenes_meanOffAncRatio <- mean(selectGenes_offAncRatio, na.rm=T)
      selectGenes_treeScore <- get_tree_score(GO_list = selectGenes_GOmin, go_graph = GO_g, log_file=logFile)
      
    } else {
      selectGenes_meanOffAncRatio <- NA
      selectGenes_treeScore <- NA
      
    }
    
    
    data.frame(
      dataset = curr_ds, 

      nSelectTADs = nSelectTADs,
      nSelectTADs_genes = nSelectTADs_genes,
      nSelectGenes = nSelectGenes,
      
      selectTADs_genes_meanOffAncRatio=selectTADs_genes_meanOffAncRatio,
      selectTADs_genes_treeScore=selectTADs_genes_treeScore,
      
      selectGenes_meanOffAncRatio=selectGenes_meanOffAncRatio,
      selectGenes_treeScore=selectGenes_treeScore,
      
      stringsAsFactors = FALSE
    )
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  
  all_ds_DT <- all_ds_DT[order(all_ds_DT$selectTADs_genes_meanOffAncRatio, 
                               all_ds_DT$selectGenes_meanOffAncRatio, 
                               decreasing=T),]
  
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}


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


all_vars_toplot <- c("meanOffAncRatio", "treeScore")
ref_var <- "selectTADs_genes"
all_cmp_var <- c("selectGenes")


var_to_plot=all_vars_toplot[1]
cmp_var=all_cmp_var[1]

for(var_to_plot in all_vars_toplot){
  for(cmp_var in all_cmp_var){
    
    var1 <- paste0(ref_var,"_",var_to_plot)
    var2 <- paste0(cmp_var,"_",var_to_plot)
    
    curr_plotDT <- plotDT_m[plotDT_m$variable %in% c(var1,var2 ),]
    
    
    plotTit <- gsub("_", " ", var_to_plot)
    mySub <- paste0("pvalSelect = ", pvalSelect, 
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

all_vars_toplot <- c("meanOffAncRatio", "treeScore")
all_ref_var <- c("aucFCC", "aucCoexprDist")
all_cmp_var <- c("selectTADs_genes","selectGenes")

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
