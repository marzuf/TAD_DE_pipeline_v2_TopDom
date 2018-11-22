
test_genes <- topTADs_genes[1:10]
test_genes <- test_genes[test_genes %in% prep_GO_Data@geneAnno$ENTREZID]
test_genes_pairs <- combn(test_genes, 2)

system.time(apply(test_genes_pairs, 2, function(x) {
  geneSim(gene1=x[1],
        gene2=x[2],
        semData=hsGO, 
        measure=semMeasure)
}))

system.time(apply(test_genes_pairs, 2, function(x) {
  TopoICSim_mz(x[1],x[2],
               GO_Data = prep_GO_Data, 
               ont="MF", organism="human", drop=NULL)
}))
