#' visual_bipartite
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_bipartite
#'
#' @import Matrix
#' @import igraph
#'
#' @param opinions number of group vertices
#' @param bipartite biparitite
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

  g <- graph_from_incidence_matrix(as.matrix(bipartite))
  V(g)$type <- bipartite_mapping(g)$type

  # Scale vertex size down for large networks
  n_nodes <- vcount(g)

  # Set colors by type
  V(g)$color <- ifelse(test = V(g)$type, yes = "yellow", no = my_palette[match(opinions, levels_vec)])
  V(g)$size <- max(10, 25 * (10 / n_nodes))
  V(g)$label <- ifelse(test = V(g)$type, yes = paste("g", 1:ncol(bipartite), sep = ""), no = NA)

  # Plot with bipartite layout
  plot(g,
       layout = layout_as_bipartite,
       vertex.color = V(g)$color,
       vertex.size = V(g)$size,
       vertex.label = V(g)$label,
       main = "Bipartite Graph")

  #Sys.sleep(0.05) # pause so you can see the update
}
