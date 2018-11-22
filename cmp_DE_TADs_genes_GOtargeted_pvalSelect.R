startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GOtargeted_pvalSelect.R\n"))

# Rscript cmp_DE_TADs_genes_GOtargeted_pvalSelect.R

# DA TADs genes vs. DE genes: C) how targeted are the results
# - number of different GO slims
# - entropy of the set of GO (kind of intra-set IC)

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ontologySimilarity, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  #GO_IC data
suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology
suppressPackageStartupMessages(library(GSEABase, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # GOCollection, goSlim
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # get_ontology

data(GO_IC)
data(go)

SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- TRUE

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
source("analysis_utils.R")

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05
ontologyType <- "MF"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  if(args[1] %in% c("BP", "MF")){
    ontologyType <- args[1]
  }
}

outFold <- file.path("CMP_TADs_GENES_GOtargeted_pvalSelect", ontologyType, pvalSelect)
if(!SSHFS) system(paste0("mkdir -p ", outFold))

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
txt <- paste0("... pvalSelect:\t", pvalSelect, "\n")
printAndLog(txt, logFile)
txt <- paste0("... ontologyType:\t", ontologyType, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:10]
# all_ds <- all_ds[1]

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
      
      ##################
      ###### # of diff. go slim - selectTADs_genes
      ##################
      # retrieve GO categories of the selectTADs_genes
      
      selectTADs_genes_GO_mapping <- bitr(selectTADs_genes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db")
      selectTADs_genes_GO_mapping <- selectTADs_genes_GO_mapping[selectTADs_genes_GO_mapping$ONTOLOGY == ontologyType,]
      selectTADs_genes_GO <- unique(selectTADs_genes_GO_mapping$GO)
      selectTADs_genes_nGO <- length(selectTADs_genes_GO)
      
      selectTADs_genes_GOmin <- minimal_set(go, selectTADs_genes_GO)
      selectTADs_genes_nGOmin <- length(selectTADs_genes_GOmin)
      
      # map the GO to their GO slim
      
      selectTADs_genes_GO_collection <- GSEABase::GOCollection(selectTADs_genes_GO)#GSEABase
      
      selectTADs_genes_GO_slim <- goSlim(selectTADs_genes_GO_collection, GOslim, ontologyType)
      
      selectTADs_genes_GOslimFreq <-  setNames(selectTADs_genes_GO_slim$Percent/100, rownames(selectTADs_genes_GO_slim))
      stopifnot(abs(1-sum(selectTADs_genes_GOslimFreq)) < 10^-6)
      txt <- paste0("... GOslim from the GO: ",  sum(selectTADs_genes_GOslimFreq > 0), "/", length(selectTADs_genes_GOslimFreq), "\n")
      printAndLog(txt, logFile)
      
      selectTADs_genes_nGOslim <- length(selectTADs_genes_GOslimFreq)
      
      # count the number of distinct go slim
      selectTADs_genes_GOslimFreq <- selectTADs_genes_GOslimFreq[selectTADs_genes_GOslimFreq > 0]
      selectTADs_genes_nGOslim <- length(selectTADs_genes_GOslimFreq)
      
      # for each of the GO slim, compute their IC
      selectTADs_genes_GOslim_IC <- -log2(selectTADs_genes_GOslimFreq)
      stopifnot(!is.na(selectTADs_genes_GOslim_IC))
      stopifnot(!is.infinite(selectTADs_genes_GOslim_IC))
      selectTADs_genes_GOslim_IC_mean <- mean(selectTADs_genes_GOslim_IC)
      stopifnot(!is.na(selectTADs_genes_GOslim_IC_mean))
      stopifnot(!is.infinite(selectTADs_genes_GOslim_IC_mean))
      
      ##################
      ###### entropy GOs - selectTADs_genes
      ##################
      # for all the GO categories from selectTADs_genes, how often they pop up
      selectTADs_genes_GOcount <- setNames(as.numeric(table(selectTADs_genes_GO_mapping$GO)),
                                           names(table(selectTADs_genes_GO_mapping$GO)))
      selectTADs_genes_GOfreq <- selectTADs_genes_GOcount/sum(selectTADs_genes_GOcount) 
      # stopifnot(sum(selectTADs_genes_GOfreq) == 1)
      selectTADs_genes_GO_IC <- -log2(selectTADs_genes_GOfreq)
      stopifnot(!is.na(selectTADs_genes_GO_IC))
      stopifnot(!is.infinite(selectTADs_genes_GO_IC))
      selectTADs_genes_GO_IC_mean <- mean(selectTADs_genes_GO_IC)
      stopifnot(!is.na(selectTADs_genes_GO_IC_mean))
      stopifnot(!is.infinite(selectTADs_genes_GO_IC_mean))
      
    } else {
      
      selectTADs_genes_nGO <- NA
      selectTADs_genes_GO_IC_mean <- NA
      
      selectTADs_genes_nGOmin <- NA
      
      selectTADs_genes_nGOslim <- NA
      selectTADs_genes_GOslim_IC_mean <- NA
    }
  
  if(nSelectGenes > 0){
    stopifnot(selectGenes %in% gene2tad_DT$entrezID)
    stopifnot(selectGenes %in% pipelineGenes)
    
    
    ##################
    ###### # of diff. go slim - selectGenes
    ##################
    # retrieve GO categories of the selectGenes
    
    selectGenes_GO_mapping <- bitr(selectGenes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db")
    
    selectGenes_GO_mapping <- selectGenes_GO_mapping[selectGenes_GO_mapping$ONTOLOGY == ontologyType,]
    selectGenes_GO <- unique(selectGenes_GO_mapping$GO)
    selectGenes_nGO <- length(selectGenes_GO)
    
    selectGenes_GOmin <- minimal_set(go, selectGenes_GO)
    selectGenes_nGOmin <- length(selectGenes_GOmin)
    
    
    # map the GO to their GO slim
    selectGenes_GO_collection <- GOCollection(selectGenes_GO)#GSEABase
    selectGenes_GO_slim <- goSlim(selectGenes_GO_collection, GOslim, ontologyType)
    selectGenes_GOslimFreq <-  setNames(selectGenes_GO_slim$Percent/100, rownames(selectGenes_GO_slim))
    stopifnot(abs(1-sum(selectGenes_GOslimFreq)) < 10^-6)
    txt <- paste0("... GOslim from the GO: ",  sum(selectGenes_GOslimFreq > 0), "/", length(selectGenes_GOslimFreq), "\n")
    printAndLog(txt, logFile)
    
    selectGenes_nGOslim <- length(selectGenes_GOslimFreq)
    
    # count the number of distinct go slim
    selectGenes_GOslimFreq <- selectGenes_GOslimFreq[selectGenes_GOslimFreq > 0]
    selectGenes_nGOslim <- length(selectGenes_GOslimFreq)
    
    # for each of the GO slim, compute their IC
    selectGenes_GOslim_IC <- -log2(selectGenes_GOslimFreq)
    stopifnot(!is.na(selectGenes_GOslim_IC))
    stopifnot(!is.infinite(selectGenes_GOslim_IC))
    selectGenes_GOslim_IC_mean <- mean(selectGenes_GOslim_IC)
    stopifnot(!is.na(selectGenes_GOslim_IC_mean))
    stopifnot(!is.infinite(selectGenes_GOslim_IC_mean))
    
    ##################
    ###### entropy GOs - selectGenes
    ##################
    # for all the GO categories from selectGenes, how often they pop up
    selectGenes_GOcount <- setNames(as.numeric(table(selectGenes_GO_mapping$GO)),
                                    names(table(selectGenes_GO_mapping$GO)))
    selectGenes_GOfreq <- selectGenes_GOcount/sum(selectGenes_GOcount) 
    # stopifnot(sum(selectGenes_GOfreq) == 1)
    selectGenes_GO_IC <- -log2(selectGenes_GOfreq)
    stopifnot(!is.na(selectGenes_GO_IC))
    stopifnot(!is.infinite(selectGenes_GO_IC))
    selectGenes_GO_IC_mean <- mean(selectGenes_GO_IC)
    stopifnot(!is.na(selectGenes_GO_IC_mean))
    stopifnot(!is.infinite(selectGenes_GO_IC_mean))
    
    
  } else {
    
    selectGenes_nGO <- NA
    selectGenes_GO_IC_mean <- NA
    
    selectGenes_nGOmin <- NA
    
    selectGenes_nGOslim <- NA
    selectGenes_GOslim_IC_mean <- NA
    
  }
    
    data.frame(
      dataset=curr_ds,
      
      nSelectTADs_genes = nSelectTADs_genes,
      
      selectTADs_genes_nGO=selectTADs_genes_nGO,
      selectTADs_genes_GO_IC_mean=selectTADs_genes_GO_IC_mean,
      
      selectTADs_genes_nGOslim=selectTADs_genes_nGOslim,
      selectTADs_genes_GOslim_IC_mean=selectTADs_genes_GOslim_IC_mean,
      
      selectTADs_genes_nGOmin=selectTADs_genes_nGOmin,
      
      nSelectGenes = nSelectGenes,
      
      selectGenes_nGO=selectGenes_nGO,
      selectGenes_GO_IC_mean=selectGenes_GO_IC_mean,
      
      selectGenes_nGOslim=selectGenes_nGOslim,
      selectGenes_GOslim_IC_mean=selectGenes_GOslim_IC_mean,
      
      selectGenes_nGOmin = selectGenes_nGOmin,
      
      stringsAsFactors = FALSE
    )
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)

  all_ds_DT <- all_ds_DT[order(all_ds_DT$selectTADs_genes_GO_IC_mean, all_ds_DT$selectGenes_GO_IC_mean, decreasing=T),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}


all_ds_DT$selectTADs_genes_nGO_ratio <- all_ds_DT$selectTADs_genes_nGO/all_ds_DT$nSelectTADs_genes
all_ds_DT$selectGenes_nGO_ratio <- all_ds_DT$selectGenes_nGO/all_ds_DT$nSelectGenes

all_ds_DT$selectTADs_genes_nGOslim_ratio <- all_ds_DT$selectTADs_genes_nGOslim/all_ds_DT$nSelectTADs_genes
all_ds_DT$selectGenes_nGOslim_ratio <- all_ds_DT$selectGenes_nGOslim/all_ds_DT$nSelectGenes

all_ds_DT$selectTADs_genes_nGOmin_ratio <- all_ds_DT$selectTADs_genes_nGOmin/all_ds_DT$nSelectTADs_genes
all_ds_DT$selectGenes_nGOmin_ratio <- all_ds_DT$selectGenes_nGOmin/all_ds_DT$nSelectGenes

#stop("-- ok --\n")


# stop("-- ok --\n")

# colnames(all_ds_DT)
# [1] "dataset"                         "selectTADs_genes_nGO"            "selectTADs_genes_GO_IC_mean"    
# [4] "selectTADs_genes_nGOslim"        "selectTADs_genes_GOslim_IC_mean" "selectGenes_nGO"                
# [7] "selectGenes_GO_IC_mean"          "selectGenes_nGOslim"             "selectGenes_GOslim_IC_mean"  


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


all_vars_toplot <- c("GO_IC_mean", "GOslim_IC_mean", "nGO", "nGOmin","nGOslim","nGO_ratio", "nGOmin_ratio","nGOslim_ratio")
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
    mySub <- paste0("ontology = ", ontologyType, 
                    " - pval. thresh = ", pvalSelect,
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
