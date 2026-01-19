#' visual_step
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_step
#'
#' @import Matrix
#' @import igraph
#'
#' @param opinions number of group vertices
#' @param RIG RIG
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step
#' # visual_step(rig_graph, opinions )
#' #---------------------------------
#' @export

visual_step <- function(RIG, opinions ){

  rig_graph <- graph_from_adjacency_matrix(as.matrix(RIG), mode = "undirected", diag = FALSE)
  V(rig_graph)$color <- ifelse(opinions == 1, "steelblue", "tomato")

  # # Scale vertex size down for large networks
  n_nodes <- vcount(rig_graph)
  # v_size <- max(3, 25 * (50 / n_nodes))  # min size = 3
  #
  # # Choose layout
  # if (n_nodes <= 20) {
  #   layout_use <- layout_in_circle(rig_graph, order = V(rig_graph))
  # } else {
  #   # Fruchterman-Reingold force-directed layout
  #   layout_use <- layout_with_fr(rig_graph, niter = 500, grid = "nogrid")
  # }


  plot(
    rig_graph,
    vertex.color = V(rig_graph)$color,
    vertex.size = max(15, 25 * (15 / n_nodes)),
    vertex.label = NA, #1:vcount(rig_graph),
    layout = layout_in_circle(rig_graph, order = V(rig_graph)),
    edge.color = "#555555",  # subtle edges
    edge.width = 3,
    edge.curved = 0.2#,  # small curvature to reduce overlap
    #main = paste("Iteration", t)
  )

  #Sys.sleep(0.05) # pause so you can see the update
}
