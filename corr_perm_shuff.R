outFold <- paste0("cmp_17f_17h_genePerm_tadShuff_2")
system(paste0("mkdir -p ", outFold))  

genePerm_dt <- eval(parse(text = load("BUILDDT_OUTPUT_FOLDER/17fpre_build_score_table_limited_pval/all_scores_DT.Rdata")))
genePerm_dt <- data.frame(genePerm_dt, stringsAsFactors = F)
rownames(genePerm_dt) <- genePerm_dt$dataset
genePerm_dt$dataset <- NULL
gpvar <- colnames(genePerm_dt)
colnames(genePerm_dt) <- paste0(colnames(genePerm_dt), "_genePerm")
genePerm_dt <- data.matrix(genePerm_dt)
stopifnot(is.numeric(genePerm_dt[1,1]))

tadShuff_dt <- eval(parse(text = load("BUILDDT_OUTPUT_FOLDER/17hpre_build_score_table_limited_pval_shuffle/all_scores_DT.Rdata")))
tadShuff_dt <- data.frame(tadShuff_dt, stringsAsFactors = F)

rownames(tadShuff_dt) <- tadShuff_dt$dataset
tadShuff_dt$dataset <- NULL
tsvar <- colnames(tadShuff_dt)
colnames(tadShuff_dt) <- paste0(colnames(tadShuff_dt), "_tadShuff")
tadShuff_dt <- data.matrix(tadShuff_dt)
stopifnot(is.numeric(tadShuff_dt[1,1]))


mergeDT <- merge(genePerm_dt, tadShuff_dt, by="row.names")

commonVar <- intersect(gpvar, tsvar)

for(i_var in commonVar) {
  
  outFile <- paste0(outFold, "/", "genePerm_tadShuff_", i_var, ".svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_genePerm")] ~  mergeDT[, paste0(i_var, "_tadShuff")],
       main = i_var,
       xlab=paste0(i_var, " - tadShuff"),
       ylab=paste0(i_var, " - genePerm"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_genePerm")],  mergeDT[, paste0(i_var, "_tadShuff")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
  outFile <- paste0(outFold, "/", "genePerm_tadShuff_", i_var, "_withLab.svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var, "_genePerm")] ~  mergeDT[, paste0(i_var, "_tadShuff")],
       main = i_var,
       xlab=paste0(i_var, " - tadShuff"),
       ylab=paste0(i_var, " - genePerm"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  text(x = mergeDT[, paste0(i_var, "_tadShuff")], 
       y = mergeDT[, paste0(i_var, "_genePerm")] ,
        labels = mergeDT$Row.names, cex=0.7)
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var, "_genePerm")],  mergeDT[, paste0(i_var, "_tadShuff")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
  
}

# TopDom - permGenes: ratioDown vs. prodSignedRatio
i_var1 <- "ratioDown_auc"
i_var2 <- "prodSignedRatio_auc"

  outFile <- paste0(outFold, "/", "genePerm_", i_var1, "_", i_var2 , "_withLab.svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var1, "_genePerm")] ~  mergeDT[, paste0(i_var2, "_genePerm")],
       main = paste0("genePerm: ", i_var1, " vs. ", i_var2),
       xlab=paste0(i_var2, " - _genePerm"),
       ylab=paste0(i_var1, " - genePerm"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  text(x = mergeDT[, paste0(i_var2, "_genePerm")], 
       y = mergeDT[, paste0(i_var1, "_genePerm")] ,
        labels = mergeDT$Row.names, cex=0.7)
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var1, "_genePerm")],  mergeDT[, paste0(i_var2, "_genePerm")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n")) 



# TopDom - tadShuff: ratioDown vs. prodSignedRatio
i_var1 <- "ratioDown_auc"
i_var2 <- "prodSignedRatio_auc"

  outFile <- paste0(outFold, "/", "tadShuff_", i_var1, "_", i_var2 , "_withLab.svg")
  svg(outFile, width=7, height=7)
  plot(mergeDT[, paste0(i_var1, "_tadShuff")] ~  mergeDT[, paste0(i_var2, "_tadShuff")],
       main = paste0("tadShuff: ", i_var1, " vs. ", i_var2),
       xlab=paste0(i_var2, " - _tadShuff"),
       ylab=paste0(i_var1, " - tadShuff"),
       cex = 0.7, pch=16)
  curve(x*1, add =T, col="grey")
  text(x = mergeDT[, paste0(i_var2, "_tadShuff")], 
       y = mergeDT[, paste0(i_var1, "_tadShuff")] ,
        labels = mergeDT$Row.names, cex=0.7)
  corTxt <- paste0("PCC = ", round(cor(mergeDT[, paste0(i_var1, "_tadShuff")],  mergeDT[, paste0(i_var2, "_tadShuff")]), 4))
  legend("topleft", legend = corTxt, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n")) 



