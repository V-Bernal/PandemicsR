#' visual_step_multi
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_step_multi
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix sparseMatrix
#' @import igraph
#'
#' @param opinions number of group vertices
#' @param RIG RIG
#' @param num_opinions number opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_multi
#' # visual_step_multi(RIG, opinions, num_opinions)
#' #---------------------------------
#' @export
visual_step_multi <- function(RIG, opinions, num_opinions){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  # Check for illegal opinions
  if(any(!(opinions %in% levels_vec))){
    stop("Unknown opinion detected: ", paste(unique(opinions), collapse=", "))
  }

  rig_graph <- graph_from_adjacency_matrix(as.matrix(RIG), mode = "undirected", diag = FALSE)

  # Map opinions to colors
  V(rig_graph)$color <- my_palette[match(opinions, levels_vec)]

  # Scale vertex size down for large networks
  n_nodes <- vcount(rig_graph)

  plot(
    rig_graph,
    vertex.color = V(rig_graph)$color,
    vertex.size = max(10, 25 * (10 / n_nodes)),
    vertex.label = NA, #1:vcount(rig_graph),
    layout = layout_in_circle(rig_graph, order = V(rig_graph)),
    edge.color = "#555555",  # subtle edges
    edge.width = 3,
    edge.curved = 0.2#,  # small curvature to reduce overlap
    #main = paste("Iteration", t)
  )

  #Sys.sleep(0.05) # pause so you can see the update
}
