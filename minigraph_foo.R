mini_graph <- inducedGraph(dag = GO_g, startNodes = selectTADs_genes_signifGOmin[selectTADs_genes_signifGOmin %in% nodes(GO_g)][1:2]) # graphNEL
plot(mini_graph)
mini_graph_rev <- reverseArch(mini_graph) # topGO
plot(mini_graph_rev)
mini_graph_rev_noAll <- removeNode("all", mini_graph_rev)
plot(mini_graph_rev_noAll)

mini_igraph <- igraph.from.graphNEL(graphNEL = mini_graph_rev_noAll)

# stopifnot(is.connected(mini_igraph))
# stopifnot(is.directed(mini_igraph))
# stopifnot(is.dag(mini_igraph))
plot(mini_igraph)
#
rootGO <- names(V(mini_igraph)[1]) # don't work with directed graph ?
#
# # Now we can find the root nodes. (either no neighbors or no incident edges)
root_idx <- which(sapply(sapply(V(mini_igraph),
                                function(x) neighbors(mini_igraph,x, mode="in")), length) == 0)
stopifnot(length(root_idx) == 1)
rootGO_bis <- names(V(mini_igraph)[root_idx])

leaves_idx <- which(sapply(sapply(V(mini_igraph),
                                  function(x) neighbors(mini_igraph, x, mode="out")), length) == 0)
leavesGO <- names(V(mini_igraph)[leaves_idx])
leavesGO
# stopifnot(rootGO == rootGO_bis)

distances(mini_igraph, v = leavesGO, to = rootGO)
nNodes <- length(V(mini_igraph))# # all is at the bottom

adjacent_vertices(mini_igraph, V(mini_igraph),mode="in")
