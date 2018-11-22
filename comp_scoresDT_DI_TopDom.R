SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

outFold <- paste0("cmp_17f_17h_DI_TopDom")
system(paste0("mkdir -p ", outFold))  

folder_TopDom <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
folder_DI <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI")

######################################################
###################################################### COMPARE GENE PERMUTATION
######################################################

genePerm_dt_TopDom <- eval(parse(text = load(paste0(folder_TopDom, "/BUILDDT_OUTPUT_FOLDER/17fpre_build_score_table_limited_pval/all_scores_DT.Rdata"))))
genePerm_dt_TopDom <- data.frame(genePerm_dt_TopDom, stringsAsFactors = F)
rownames(genePerm_dt_TopDom) <- genePerm_dt_TopDom$dataset
genePerm_dt_TopDom$dataset <- NULL
gpvar_TopDom <- colnames(genePerm_dt_TopDom)
colnames(genePerm_dt_TopDom) <- paste0(colnames(genePerm_dt_TopDom), "_genePermTopDom")
genePerm_dt_TopDom <- data.matrix(genePerm_dt_TopDom)
stopifnot(is.numeric(genePerm_dt_TopDom[1,1]))


genePerm_dt_DI <- eval(parse(text = load(paste0(folder_DI, "/BUILDDT_OUTPUT_FOLDER/17fpre_build_score_table_limited_pval/all_scores_DT.Rdata"))))
genePerm_dt_DI <- data.frame(genePerm_dt_DI, stringsAsFactors = F)
rownames(genePerm_dt_DI) <- genePerm_dt_DI$dataset
genePerm_dt_DI$dataset <- NULL
gpvar_DI <- colnames(genePerm_dt_DI)
colnames(genePerm_dt_DI) <- paste0(colnames(genePerm_dt_DI), "_genePermDI")
genePerm_dt_DI <- data.matrix(genePerm_dt_DI)
stopifnot(is.numeric(genePerm_dt_DI[1,1]))



mergeDT <- merge(genePerm_dt_TopDom, genePerm_dt_DI, by="row.names")

commonVar <- intersect(gpvar_TopDom, gpvar_DI)
stopifnot(length(commonVar) > 0)

for(i_var in commonVar) {
  outFile <- paste0(outFold, "/", "genePerm_TopDom_DI_", i_var, ".svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_genePermTopDom")] ~  mergeDT[, paste0(i_var, "_genePermDI")],
       main = paste0(i_var, " -  gene permuation"),
       xlab=paste0(i_var, " - DI"),
       ylab=paste0(i_var, " - TopDom"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_genePermTopDom")] ,  mergeDT[, paste0(i_var, "_genePermDI")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
  outFile <- paste0(outFold, "/", "genePerm_TopDom_DI_", i_var, "_withLab.svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_genePermTopDom")] ~  mergeDT[, paste0(i_var, "_genePermDI")],
       main = paste0(i_var, " -  gene permuation"),
       xlab=paste0(i_var, " - DI"),
       ylab=paste0(i_var, " - TopDom"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  text(x = mergeDT[, paste0(i_var, "_genePermDI")], 
       y = mergeDT[, paste0(i_var, "_genePermTopDom")] ,
       labels = mergeDT$Row.names, cex=0.7)
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_genePermTopDom")] , mergeDT[, paste0(i_var, "_genePermDI")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
}


######################################################
###################################################### COMPARE TAD SHUFFLING
######################################################

tadShuff_dt_TopDom <- eval(parse(text = load(paste0(folder_TopDom, "/BUILDDT_OUTPUT_FOLDER/17hpre_build_score_table_limited_pval_shuffle/all_scores_DT.Rdata"))))
tadShuff_dt_TopDom <- data.frame(tadShuff_dt_TopDom, stringsAsFactors = F)

rownames(tadShuff_dt_TopDom) <- tadShuff_dt_TopDom$dataset
tadShuff_dt_TopDom$dataset <- NULL
tsvar <- colnames(tadShuff_dt_TopDom)
colnames(tadShuff_dt_TopDom) <- paste0(colnames(tadShuff_dt_TopDom), "_tadShuffTopDom")
tadShuff_dt_TopDom <- data.matrix(tadShuff_dt_TopDom)
stopifnot(is.numeric(tadShuff_dt_TopDom[1,1]))



tadShuff_dt_DI <- eval(parse(text = load(paste0(folder_DI, "/BUILDDT_OUTPUT_FOLDER/17hpre_build_score_table_limited_pval_shuffle/all_scores_DT.Rdata"))))
tadShuff_dt_DI <- data.frame(tadShuff_dt_DI, stringsAsFactors = F)

rownames(tadShuff_dt_DI) <- tadShuff_dt_DI$dataset
tadShuff_dt_DI$dataset <- NULL
tsvar <- colnames(tadShuff_dt_DI)
colnames(tadShuff_dt_DI) <- paste0(colnames(tadShuff_dt_DI), "_tadShuffDI")
tadShuff_dt_DI <- data.matrix(tadShuff_dt_DI)
stopifnot(is.numeric(tadShuff_dt_DI[1,1]))

mergeDT <- merge(tadShuff_dt_TopDom, tadShuff_dt_DI, by="row.names")


for(i_var in commonVar) {
  outFile <- paste0(outFold, "/", "tadShuff_TopDom_DI_", i_var, ".svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_tadShuffTopDom")] ~  mergeDT[, paste0(i_var, "_tadShuffDI")],
       main = paste0(i_var, " -  TAD shuffling"),
       xlab=paste0(i_var, " - DI"),
       ylab=paste0(i_var, " - TopDom"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_tadShuffTopDom")] ,  mergeDT[, paste0(i_var, "_tadShuffDI")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
  outFile <- paste0(outFold, "/", "tadShuff_TopDom_DI_", i_var, "_withLab.svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_tadShuffTopDom")] ~  mergeDT[, paste0(i_var, "_tadShuffDI")],
       main = paste0(i_var, " -  TAD shuffling"),
       xlab=paste0(i_var, " - DI"),
       ylab=paste0(i_var, " - TopDom"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  text(x = mergeDT[, paste0(i_var, "_tadShuffDI")], 
       y = mergeDT[, paste0(i_var, "_tadShuffTopDom")] ,
       labels = mergeDT$Row.names, cex=0.7)
  
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_tadShuffTopDom")] ,  mergeDT[, paste0(i_var, "_tadShuffDI")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
}



