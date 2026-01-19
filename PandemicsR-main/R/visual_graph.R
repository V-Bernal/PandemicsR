#' visual_graph
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_graph
#'
#' @param RIG RIG
#' @param num_opinions num_opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_graph
#' # visual_graph(RIG, num_opinions)
#' #---------------------------------
#' @export
visual_graph <- function(RIG, num_opinions){

  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  rig_graph <- graph_from_adjacency_matrix(as.matrix(RIG), mode = "undirected", diag = FALSE)

  par(mfrow = c(4, 1))

  deg <- degree(rig_graph)
  hist(deg, breaks = "FD",
       main = "Degree distribution",
       xlab = "Degree")

  plot(sort(deg, decreasing = TRUE),
       log = "xy",
       main = "Degree rank plot",
       xlab = "Rank", ylab = "Degree")

  transitivity(rig_graph, type = "global")
  hist(transitivity(rig_graph, type = "local", isolates = "zero"),
       main = "Local clustering coefficient",
       xlab = "C")

  comp <- components(rig_graph)
  barplot(comp$csize,
          main = "Component size distribution",
          xlab = "Component", ylab = "Size")

}
