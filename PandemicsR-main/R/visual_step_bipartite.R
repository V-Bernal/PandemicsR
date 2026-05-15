#' visual_bipartite
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_bipartite
#'
#' @param opinions number of group vertices
#' @param bipartite bipartite
#' @param num_opinions number opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_bipartite
#' # visual_bipartite(bipartite, opinions, num_opinions)
#' #---------------------------------
#' @export

visual_bipartite <- function(bipartite, opinions, num_opinions){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  # Check for illegal opinions
  if(any(!(opinions %in% levels_vec))){
    stop("Unknown opinion detected: ", paste(unique(opinions), collapse=", "))
  }

  bipartite_graph <- Matrix::Matrix(bipartite != 0, sparse = TRUE)
  n_individuals <- nrow(bipartite_graph)
  n_groups <- ncol(bipartite_graph)
  edge_count <- Matrix::nnzero(bipartite_graph)

  if (n_individuals > 120 || n_groups > 20 || edge_count > 1500) {
    x_ind <- seq(0.02, 0.98, length.out = n_individuals)
    x_grp <- seq(0.05, 0.95, length.out = max(1L, n_groups))
    edges <- Matrix::summary(bipartite_graph)

    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    if (nrow(edges) > 0) {
      segments(
        x0 = x_ind[edges$i],
        y0 = 0.72,
        x1 = x_grp[edges$j],
        y1 = 0.28,
        col = grDevices::adjustcolor("gray45", alpha.f = 0.18),
        lwd = 0.5
      )
    }
    points(
      x_ind,
      rep(0.72, n_individuals),
      pch = 16,
      col = my_palette[match(opinions, levels_vec)],
      cex = max(0.3, 8 / sqrt(n_individuals))
    )
    points(
      x_grp,
      rep(0.28, n_groups),
      pch = 21,
      bg = "yellow",
      col = "goldenrod4",
      cex = max(0.8, 4 / sqrt(max(1, n_groups)))
    )
    if (n_groups <= 30) {
      text(x_grp, rep(0.23, n_groups), labels = paste0("g", seq_len(n_groups)), cex = 0.7)
    }
    box()
    title(main = "Bipartite Graph", sub = "Simplified view for larger networks")
    return(invisible(NULL))
  }

  g <- igraph::graph_from_biadjacency_matrix(bipartite_graph)

  # Scale vertex size down for large networks
  n_nodes <- igraph::vcount(g)

  # Set colors by type
  vertex_type <- igraph::V(g)$type
  vertex_color <- ifelse(test = vertex_type, yes = "yellow", no = my_palette[match(opinions, levels_vec)])
  vertex_size <- max(10, 25 * (10 / n_nodes))
  vertex_label <- ifelse(test = vertex_type, yes = paste("g", 1:ncol(bipartite), sep = ""), no = NA)

  # Plot with bipartite layout
  plot(g,
       layout = igraph::layout_as_bipartite,
       vertex.color = vertex_color,
       vertex.size = vertex_size,
       vertex.label = vertex_label,
       edge.color = grDevices::adjustcolor("gray45", alpha.f = 0.45),
       main = "Bipartite Graph")

  #Sys.sleep(0.05) # pause so you can see the update
}
