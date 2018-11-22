startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_summaryDT.R\n"))

suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# Rscript cmp_DE_TADs_genes_summaryDT.R

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2,40))

outFold <- "CMP_DE_TADs_SUMMARY"
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 500, 8)
myWidth <- ifelse(plotType == "png", 600, 10)

caller <- "TopDom"
pipOutFold <- "OUTPUT_FOLDER"
script170_name <- "170_score_auc_pval_withShuffle"

barcolor <- "navy"

source("analysis_utils.R")

# retrieve data generated from the following scripts:
# - cmp_DE_TADs_genes_GO.R
#     -> /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_TADs_GENES_GO/all_ds_DT.Rdata
# [1] "dataset"                                   
# [2] "topTADs_genes_topGo"                       
# [3] "topGenes_manyAsTopPercent_topGo"           
# [4] "topGenes_manyAsTopTADs_topGo"              
# [5] "nTopPercentTADs"                           
# [6] "nTopPercentTADs_genes"                     
# [7] "nTopPercentGenes"                          
# [8] "intersectNbrGO_topTADs_topGenesManyAsTopPercent"
# [9] "intersectNbrGO_topTADs_topGenesManyAsTopTADs"   
goFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_TADs_GENES_GO")
goFile <- file.path(goFold, "all_ds_DT.Rdata")
stopifnot(file.exists(goFile))

# - cmp_DE_TADs_genes_intersectSignif.R
#     -> /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_TADs_GENES_INTERSECT/all_ds_DT.Rdata
# [1] "dataset"               
# [2] "nTopPercentTADs"       
# [3] "nTopPercentTADs_genes"
# [4] "nTopPercentGenes"      
# [5] "intersectNbrGenes_nbrPercent"   
# [6] "intersectNbrGenes_sameNbr"     
intersectFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_TADs_GENES_INTERSECT")
intersectFile <- file.path(intersectFold, "all_ds_DT.Rdata")
stopifnot(file.exists(intersectFile))

# - cmp_DE_TADs_genes_sum_ranks.R
#     -> /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_SUM_RANKS/all_ds_all_TADs_DT.Rdata
# [1] "dataset"        
# [2] "TAD"            
# [3] "TAD_rank"       
# [4] "TAD_genes"     
# [5] "TAD_genes_rank" 
# [6] "nGenes"         
# [7] "sumGeneRank"    
# [8] "genes_meanRank"   
rankFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CMP_SUM_RANKS")
rankFile <- file.path(rankFold, "all_ds_all_TADs_DT.Rdata")
stopifnot(file.exists(rankFile))

go_DT <- eval(parse(text = load(goFile)))
stopifnot(colnames(go_DT) == c("dataset","topTADs_genes_topGo","topGenes_manyAsTopPercent_topGo",
                               "topGenes_manyAsTopTADs_topGo","nTopPercentTADs","nTopPercentTADs_genes",
                               "nTopPercentGenes","nbrGO_topTADs_genes",
                               "nbrGO_topGenes_manyAsTopPercent","nbrGO_topGenes_manyAsTopTADs",
                               "intersectNbrGO_topTADs_topGenesManyAsTopPercent","intersectNbrGO_topTADs_topGenesManyAsTopTADs"))
go_DT$dataset <- as.character(go_DT$dataset)

intersect_DT <- eval(parse(text = load(intersectFile)))
stopifnot(colnames(intersect_DT) == c("dataset","nTopPercentTADs","nTopPercentTADs_genes","nTopPercentGenes","intersectNbrGenes_nbrPercent","intersectNbrGenes_sameNbr"))
intersect_DT$dataset <- as.character(intersect_DT$dataset)

rank_DT <- eval(parse(text = load(rankFile)))
stopifnot(colnames(rank_DT) == c("dataset", "TAD", "TAD_rank","TAD_genes", "TAD_genes_rank","nGenes","sumGeneRank","genes_meanRank"))
rank_DT$dataset <- as.character(rank_DT$dataset)
rank_DT$diffRank <- abs(rank_DT$genes_meanRank - rank_DT$TAD_rank)
# aggregate to have 1 value per dataset
agg_rank_DT <- aggregate(diffRank ~ dataset, data = rank_DT, FUN=mean)
agg_rank_DT$dataset <- as.character(agg_rank_DT$dataset)

stopifnot(!length(setdiff(intersect_DT$dataset,go_DT$dataset)))
stopifnot(!length(setdiff(intersect_DT$dataset,rank_DT$dataset)))
stopifnot(!length(setdiff(intersect_DT$dataset,agg_rank_DT$dataset)))
stopifnot(!length(setdiff(go_DT$dataset,rank_DT$dataset)))
stopifnot(!length(setdiff(go_DT$dataset,agg_rank_DT$dataset)))

stopifnot(!duplicated(intersect_DT$dataset))
stopifnot(!duplicated(go_DT$dataset))
stopifnot(!duplicated(agg_rank_DT$dataset))

# build common table
go_intersect_DT <- merge(go_DT, intersect_DT, by=c("dataset", "nTopPercentTADs","nTopPercentTADs_genes", "nTopPercentGenes"))
stopifnot(!duplicated(go_intersect_DT$dataset))
stopifnot(nrow(go_intersect_DT) == nrow(intersect_DT))
stopifnot(nrow(go_intersect_DT) == nrow(go_DT))
go_intersect_DT <- go_intersect_DT[, c("dataset", 
                                       "intersectNbrGO_topTADs_topGenesManyAsTopPercent", 
                                       "intersectNbrGO_topTADs_topGenesManyAsTopTADs", 
                                       "intersectNbrGenes_nbrPercent", 
                                       "intersectNbrGenes_sameNbr")]

summary_DT <- merge(go_intersect_DT, agg_rank_DT, by ="dataset")
summary_DT$dataset <- as.character(summary_DT$dataset)
stopifnot(!duplicated(summary_DT$dataset))
stopifnot(nrow(summary_DT) == nrow(agg_rank_DT))
stopifnot(nrow(summary_DT) == nrow(go_intersect_DT))

outFile <- file.path(outFold, "summary_DT.Rdata")
save(summary_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

var_to_plot <- colnames(summary_DT)[colnames(summary_DT) != "dataset"]

all_datasets <- summary_DT$dataset
# retrieved the AUC ratio for all the dataset
all_auc <- foreach(curr_dataset = all_datasets) %dopar% {
  
  aucFCC_file <- file.path(pipOutFold, curr_dataset, script170_name, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_dataset, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  
  list(aucFCC = aucFCC, aucCoexprDist = aucCoexprDist)
  
}
names(all_auc) <- all_datasets


curr_var <- "intersectNbrGO_topTADs_topGenesManyAsTopTADs"

var_to_plot_title <- c(
"intersectNbrGO_topTADs_topGenesManyAsTopPercent" = "# GO intersect (# genes in top%)",
"intersectNbrGO_topTADs_topGenesManyAsTopTADs" = "# GO intersect (# genes = # topTADs genes)",
"intersectNbrGenes_nbrPercent" = "# genes intersect (# genes in top%)",
"intersectNbrGenes_sameNbr" = "# genes intersect (# genes = # topTADs genes)",
"diffRank" = "TAD-genes rank difference"
)

stopifnot(var_to_plot %in% names(var_to_plot_title))

for(curr_var in var_to_plot) {
  
  plotDT <- summary_DT[,c("dataset", curr_var)]
  plotDT <- plotDT[order(plotDT[,curr_var], plotDT[,"dataset"], decreasing = T),]
  plotDT$dataset <- factor(plotDT$dataset, levels = as.character(plotDT$dataset))
  
  
  p_var <- ggplot(plotDT, aes_string(x = "dataset", y = curr_var)) +
    geom_bar(stat="identity", fill = barcolor)+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(var_to_plot_title[curr_var]),
                       breaks = scales::pretty_breaks(n = 5))+
    theme( # Increase size of axis lines
      # strip.text = element_text(size = 12),
      # top, right, bottom and left
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=10),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
    )
  
  outFile <- file.path(outFold, paste0("barplot_datasets_", curr_var, ".", plotType))
  ggsave(p_var, filename = outFile, height = myHeight, width = myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  plotDT$aucFCC <- sapply(plotDT$dataset, function(x) all_auc[[x]][["aucFCC"]])
  plotDT$aucCoexprDist <- sapply(plotDT$dataset, function(x) all_auc[[x]][["aucCoexprDist"]])
  
  curr_var2 <- "aucFCC"
  for(curr_var2 in c("aucFCC", "aucCoexprDist")) {

    myx <- plotDT[,curr_var2]
    myy <- plotDT[,curr_var]

    outFile <- file.path(outFold, paste0("scatter_", curr_var, "_vs_", curr_var2, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(x = myx,
         y = myy,
         ylab = paste0(curr_var),
         xlab = paste0(curr_var2),
         pch=16, cex = 0.7,
         main=paste0(var_to_plot_title[curr_var], " vs. ", curr_var2))
    addCorr(x=myx, y=myy, legPos="topleft", bty="n")
    add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
    foo <- dev.off()    
    cat(paste0("... written: ", outFile, "\n"))


    outFile <- file.path(outFold, paste0("scatter_", curr_var, "_vs_", curr_var2, "_withLabs.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(x = myx,
         y = myy,
         ylab = paste0(curr_var),
         xlab = paste0(curr_var2),
         pch=16, cex = 0.7,
         main=paste0(var_to_plot_title[curr_var], " vs. ", curr_var2))
    addCorr(x=myx, y=myy, legPos="topleft", bty="n")
    add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
    text(x=myx, y=myy, labels=plotDT[,"dataset"], cex=0.7)
    foo <- dev.off()    
    cat(paste0("... written: ", outFile, "\n"))




  }
  
  
}



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
