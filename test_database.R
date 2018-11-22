suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

registerDoMC(30)

my_genes <- c("835",  "5261" ,"241" )

combineSemSimMethod = "BMA"

semSimMetric = "Wang"
  
hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)
  
cat("... test in lapply\n")
topTADs_semSim <- lapply(c(1:100), function(x) {
  # cat("... compute TRUE semantic similarity for TAD:", x, "\n")
  # tad_genes <- topTADs_genes[[x]]
  tad_semSim <-  mgeneSim(genes=my_genes,
                          semData=hsGO, 
                          combine=combineSemSimMethod,
                          measure=semSimMetric,
                          verbose=FALSE)
})

cat("***lapply done\n")

cat("... test in foreach\n")

topTADs_semSim <- foreach(i=c(1:100)) %dopar% {
  # cat("... compute TRUE semantic similarity for TAD:", x, "\n")
  # tad_genes <- topTADs_genes[[x]]
  tad_semSim <-  mgeneSim(genes=my_genes,
                          semData=hsGO, 
                          combine=combineSemSimMethod,
                          measure=semSimMetric,
                          verbose=FALSE)
}


cat("***foreach done\n")

cat("... test foreach in lapply \n")

topTADs_semSim <- lapply(c(1:100), function(x) {
  # cat("... compute TRUE semantic similarity for TAD:", x, "\n")
  # tad_genes <- topTADs_genes[[x]]
  tad_semSim <-  mgeneSim(genes=my_genes,
                          semData=hsGO, 
                          combine=combineSemSimMethod,
                          measure=semSimMetric,
                          verbose=FALSE)
  
  topTADs_semSim <- foreach(i=c(1:100)) %dopar% {
    # cat("... compute TRUE semantic similarity for TAD:", x, "\n")
    # tad_genes <- topTADs_genes[[x]]
    tad_semSim <-  mgeneSim(genes=my_genes,
                            semData=hsGO, 
                            combine=combineSemSimMethod,
                            measure=semSimMetric,
                            verbose=FALSE)
  }
})

