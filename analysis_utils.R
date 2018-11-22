add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}

addCorr <- function(x, y, legPos="topright", corMet="pearson", ...) {
  corMet <- tolower(corMet)
  stopifnot(corMet %in% c("pearson", "kendall", "spearman"))
  ct <- cor.test(x,y, method = corMet)
  corCoeff <- ct$estimate
  corPval <- ct$p.value
  legTxt <- paste0(paste0(toupper(substr(corMet,1,1)), "CC"), " = ", round(corCoeff, 4), "\n", "(p-val = ", sprintf("%2.2e", corPval), ")")
  legend(legPos, legend = legTxt, ...)
}

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}


