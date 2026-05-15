#' visual_step_multi
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_step_multi
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

  n_nodes <- length(opinions)
  edge_matrix <- Matrix::triu(RIG != 0, k = 1)
  edge_count <- Matrix::nnzero(edge_matrix)

  if (n_nodes > 250 || edge_count > 5000) {
    theta <- seq(0, 2 * pi, length.out = n_nodes + 1L)[-(n_nodes + 1L)]
    x <- cos(theta)
    y <- sin(theta)
    edges <- Matrix::summary(edge_matrix)

    plot.new()
    plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
    if (nrow(edges) > 0) {
      segments(
        x0 = x[edges$i],
        y0 = y[edges$i],
        x1 = x[edges$j],
        y1 = y[edges$j],
        col = grDevices::adjustcolor("#555555", alpha.f = 0.12),
        lwd = 0.6
      )
    }
    points(
      x,
      y,
      pch = 16,
      col = my_palette[match(opinions, levels_vec)],
      cex = max(0.35, 8 / sqrt(n_nodes))
    )
    box()
    title(main = "RIG", sub = "Simplified view for larger networks")
    return(invisible(NULL))
  }

  rig_graph <- igraph::graph_from_adjacency_matrix(RIG, mode = "undirected", diag = FALSE)

  # Map opinions to colors
  vertex_color <- my_palette[match(opinions, levels_vec)]

  # Scale vertex size down for large networks
  n_nodes <- igraph::vcount(rig_graph)

  plot(
    rig_graph,
    vertex.color = vertex_color,
    vertex.size = max(10, 25 * (10 / n_nodes)),
    vertex.label = NA, #1:vcount(rig_graph),
    layout = igraph::layout_in_circle(rig_graph, order = igraph::V(rig_graph)),
    edge.color = "#555555",  # subtle edges
    edge.width = max(0.8, min(3, 120 / max(1, n_nodes))),
    edge.curved = 0#,  # straight edges are faster and cleaner here
    #main = paste("Iteration", t)
  )

  #Sys.sleep(0.05) # pause so you can see the update
}
