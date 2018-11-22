#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================

get_offAncRatio <- function(GO_list, go_data_off, go_data_anc) {
  
  # offspring-ancestor ratio as a measure of specificity
  genes_offAncRatio <- sapply(GO_list, function(x) {
    nOffspring <- length(go_data_off[[x]])
    nAncestor <- length(go_data_anc[[x]])
    if(nOffspring > 0 & nAncestor > 0) {
      return(nOffspring/(nOffspring + nAncestor))
    }else{
      return(NA)
    }
  })
  return(genes_offAncRatio)
}

#    topTADs_genes_GO_mapping <- bitr(topTADs_genes, fromType="ENTREZID", toType=c("GO"), OrgDb="org.Hs.eg.db") # clusterProfiler
#    topTADs_genes_GO_mapping <- topTADs_genes_GO_mapping[topTADs_genes_GO_mapping$ONTOLOGY == ontologyType,]
#    topTADs_genes_GO <- unique(topTADs_genes_GO_mapping$GO)
#    topTADs_genes_GOmin <- minimal_set(go, topTADs_genes_GO) # ontologyIndex
#    
#    # offspring-ancestor ratio as a measure of specificity
#    topTADs_genes_offAncRatio <- sapply(topTADs_genes_GOmin, function(x) {
#      nOffspring <- length(GO_OFFSPRING[[x]])
#      nAncestor <- length(GO_ANCESTOR[[x]])
#      if(nOffspring > 0 & nAncestor > 0) {
#        return(nOffspring/(nOffspring + nAncestor))
#      }else{
#        return(NA)
#      }
#    })

#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================

get_tree_score <- function(GO_list, go_graph, log_file=""){
 
 stopifnot(any(GO_list %in% nodes(go_graph)))
 txt <- paste0("...... matching GOmin in graph nodes: ", sum(GO_list %in% nodes(go_graph)), "/", length(GO_list), "\n")
 printAndLog(txt, log_file)
 
 # retrieve the full graph for the GO
 GO_list_graph_tmp <- inducedGraph(dag = go_graph, startNodes = GO_list[GO_list %in% nodes(go_graph)]) # graphNEL
 # !!! NEED TO REVERSE THE GRAPH !!!
 GO_list_graph_tmp_rev <- reverseArch(GO_list_graph_tmp) # topGO
 GO_list_graph_tmp_rev_noAll <- removeNode("all", GO_list_graph_tmp_rev)
 
 GO_list_graph <- igraph.from.graphNEL(graphNEL = GO_list_graph_tmp_rev_noAll)
 
 stopifnot(is.connected(GO_list_graph))
 stopifnot(is.directed(GO_list_graph))
 stopifnot(is.dag(GO_list_graph))
 
 nNodes <- length(V(GO_list_graph))
 
 # Now we can find the root nodes. (either no neighbors or no incident edges)
 root_idx <- which(sapply(sapply(V(GO_list_graph), 
                                 function(x) neighbors(GO_list_graph, x, mode="in")), length) == 0)
 stopifnot(length(root_idx) == 1)
 rootGO <- names(V(GO_list_graph)[root_idx])
 txt <- paste0("...... found root GO: ", rootGO, "\n")
 printAndLog(txt, log_file)
 
 # distances from all leaves to the root
 leaves_idx <- which(sapply(sapply(V(GO_list_graph), 
                                   function(x) neighbors(GO_list_graph, x, mode="out")), length) == 0)
 txt <- paste0("... found # leaves GO: ", length(leaves_idx),  "/", length(GO_list), "\n")
 printAndLog(txt, log_file)
 
 leavesGO <- names(V(GO_list_graph)[leaves_idx])
 nLeaves <- length(leavesGO)
 
 all_leaves_dist <- distances(GO_list_graph, v = leavesGO, to = rootGO)
 stopifnot(ncol(all_leaves_dist) == 1)
 stopifnot(nrow(all_leaves_dist) == nLeaves)
 stopifnot(colnames(all_leaves_dist) == rootGO)
 
 meanLeavesDist <- mean(all_leaves_dist[,1])
 stopifnot(!is.na(meanLeavesDist))
 
 tree_score <- 1 - meanLeavesDist/(nNodes -1)
 
 stopifnot(tree_score >= 0)
 stopifnot(tree_score <= 1)
 
 return(tree_score)
}



#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================






my_plot_GO_graph <- function(GOgraph, GOids, graphTit=""){
  
  par(mfrow=c(1,2))
  
  gNEL <- inducedGraph(dag = GOgraph, startNodes = GOids)
  # !!! need to reverse the graph !!!
  gNEL_rev <- reverseArch(gNEL) # topGO
  gNEL_rev_noAll <- removeNode("all", gNEL_rev) # not needed ? -> do not want this root
  ig <- igraph.from.graphNEL(graphNEL = gNEL_rev_noAll)
  plot(gNEL_rev_noAll)
  
  ig_nNodes <- length(V(ig))
  
  # density
  ig_density <- edge_density(ig)
  
  # edge connectivity (adhesion)
  # ig_undir <- as.undirected(ig)
  ig_adhesion <- edge_connectivity(ig)
  
  # vertex connectivity (cohesion)
  ig_cohesion <- vertex_connectivity(ig)
  # diameter
  ig_diameter <- diameter(ig)
  
  # mean geodesic distance
  ig_meanDist <- mean_distance(ig, directed = FALSE) # directed TRUE or FALSE ???
  ig_meanDistDir <- mean_distance(ig, directed = TRUE) # directed TRUE or FALSE ???
  
  plot(ig, layout=layout_as_tree)
  title(main="selectTADs_genes subgraph",
        sub = paste0(
          "nNodes = ", ig_nNodes, ", \t",
          "adhesion = ", ig_adhesion, ", \t",
          "cohesion = ",ig_cohesion , ", \t",
          "diameter = ", ig_diameter, ", \n",
          "density = ", round(ig_density, 4), ", \t",
          "meanDist = ", round(ig_meanDist, 4) , ", \t",
          "meanDistDir = ", round(ig_meanDistDir, 4), "\t"
        ))
  mtext(paste0(GOids), side = 3, line = -2, font = 3)
  
  
}





#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================

#  # tree score as a measure of targeted    
#    # GO_g defined outside the loop as 
#    # GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)
#    stopifnot(any(topTADs_genes_GOmin %in% nodes(GO_g)))
#    txt <- paste0("...... matching GOmin in graph nodes: ", sum(topTADs_genes_GOmin %in% nodes(GO_g)), "/", length(topTADs_genes_GOmin), "\n")
#    printAndLog(txt, logFile)
#    
#    # retrieve the full graph for the GO
#    topTADs_genes_GOmin_graph_tmp <- inducedGraph(dag = GO_g, startNodes = topTADs_genes_GOmin[topTADs_genes_GOmin %in% nodes(GO_g)]) # graphNEL
#    # !!! NEED TO REVERSE THE GRAPH !!!
#    topTADs_genes_GOmin_graph_tmp_rev <- reverseArch(topTADs_genes_GOmin_graph_tmp) # topGO
#    topTADs_genes_GOmin_graph_tmp_rev_noAll <- removeNode("all", topTADs_genes_GOmin_graph_tmp_rev)
#    
#    topTADs_genes_GOmin_graph <- igraph.from.graphNEL(graphNEL = topTADs_genes_GOmin_graph_tmp_rev_noAll)
#    
#    stopifnot(is.connected(topTADs_genes_GOmin_graph))
#    stopifnot(is.directed(topTADs_genes_GOmin_graph))
#    stopifnot(is.dag(topTADs_genes_GOmin_graph))
#    
#    # Now we can find the root nodes. (either no neighbors or no incident edges)
#    root_idx <- which(sapply(sapply(V(topTADs_genes_GOmin_graph), 
#                                    function(x) neighbors(topTADs_genes_GOmin_graph, x, mode="in")), length) == 0)
#    stopifnot(length(root_idx) == 1)
#    rootGO <- names(V(topTADs_genes_GOmin_graph)[root_idx])
#    
#    txt <- paste0("...... found root GO: ", rootGO, "\n")
#    printAndLog(txt, logFile)
#    
#    
#    # distances from all leaves to the root
#    leaves_idx <- which(sapply(sapply(V(topTADs_genes_GOmin_graph), 
#                                    function(x) neighbors(topTADs_genes_GOmin_graph, x, mode="out")), length) == 0)
#    txt <- paste0("... found # leaves GO: ", length(leaves_idx),  "/", length(topTADs_genes_GOmin), "\n")
#    printAndLog(txt, logFile)
#    leavesGO <- names(V(topTADs_genes_GOmin_graph)[leaves_idx])
#    
#    nLeaves <- length(leavesGO)
#    
#    all_leaves_dist <- distances(topTADs_genes_GOmin_graph, v = leavesGO, to = rootGO)
#    stopifnot(ncol(all_leaves_dist) == 1)
#    stopifnot(nrow(all_leaves_dist) == nLeaves)
#    stopifnot(colnames(all_leaves_dist) == rootGO)
#    
#    meanLeavesDist <- mean(all_leaves_dist[,1])
#    stopifnot(!is.na(meanLeavesDist))
#    
#    nNodes <- length(V(topTADs_genes_GOmin_graph))
#    
#    tree_score <- 1 - meanLeavesDist/(nNodes -1)
    
    
    #######################################################
    ####################################################### FOO SMALL GRAPH
    #######################################################
    # # # investigate how looks like for small graph 
        # mini_graph <- inducedGraph(dag = GO_g, startNodes = topTADs_genes_GOmin[topTADs_genes_GOmin %in% nodes(GO_g)][1:2]) # graphNEL
        # plot(mini_graph)
        # # # all is at the bottom
    # mini_graph_rev <- reverseArch(mini_graph) # topGO
    # plot(mini_graph_rev)
    # mini_graph_rev_noAll <- removeNode("all", mini_graph_rev)
    # plot(mini_graph_rev_noAll)
    # 
    # mini_igraph <- igraph.from.graphNEL(graphNEL = mini_graph_rev_noAll)
    # 
    # # stopifnot(is.connected(mini_igraph))
    # # stopifnot(is.directed(mini_igraph))
    # # stopifnot(is.dag(mini_igraph))
    # plot(mini_igraph)
    # #     
    # rootGO <- names(V(mini_igraph)[1]) # don't work with directed graph ?
    # # 
    # # # Now we can find the root nodes. (either no neighbors or no incident edges)
    # root_idx <- which(sapply(sapply(V(mini_igraph),
    #                     function(x) neighbors(mini_igraph,x, mode="in")), length) == 0)
    # stopifnot(length(root_idx) == 1)
    # rootGO_bis <- names(V(mini_igraph)[root_idx])
    # 
    # leaves_idx <- which(sapply(sapply(V(mini_igraph), 
    #                                   function(x) neighbors(mini_igraph, x, mode="out")), length) == 0)
    # leavesGO <- names(V(mini_igraph)[leaves_idx])
    # leavesGO
    # # stopifnot(rootGO == rootGO_bis)
    # 
    # distances(mini_igraph, v = leavesGO, to = rootGO)
    # nNodes <- length(V(mini_igraph))
    ##############################################################################################################
    


