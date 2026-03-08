#' visual_histo_pergroup
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_histo_pergroup
#'
#' @import Matrix
#' @import igraph
#'
#' @references \url{}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @param opinions opinions
#' @param num_opinions number of opinions
#' @param bipartite bipartite
#'
#' @return visual_histo_pergroup \code{x}.
#'
#' @examples
#' #---------------------------------
#' # visual_histo_pergroup
#' # visual_histo_pergroup_multi(opinions, num_opinions, bipartite)
#' #---------------------------------
#' @export
visual_histo_pergroup <- function(opinions, num_opinions, members) {
  
  # Helper objects
  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)
  
  m <- length(members)
  
  # Build frequency matrix directly
  freq_mat <- sapply(seq_len(m), function(g) {
    
    ids <- members[[g]]
    
    if (length(ids) == 0) {
      return(rep(0, length(levels_vec)))
    }
    
    tab <- table(factor(opinions[ids], levels = levels_vec))
    as.numeric(tab)
  })
  
  rownames(freq_mat) <- levels_vec
  colnames(freq_mat) <- paste("group", seq_len(m))
  
  # Plot
  barplot(freq_mat,
          names.arg = colnames(freq_mat),
          col = my_palette,
          main = "Opinions",
          beside = FALSE)
}



# visual_histo_pergroup <- function(opinions, num_opinions, bipartite) {
# 
#   # Helper objects
#   levels_vec <- get_levels_vec(num_opinions)   # vector of opinion levels
#   my_palette <- get_palette(num_opinions)      # colors for levels
# 
#   # Multiply matrices to get weighted opinions
#   #tb0 <- as.matrix(bipartite0 * opinion_history[, 1])
#   #tb <- as.matrix(bipartite * opinions)
#   tb <- as.matrix(sweep(bipartite, 1, opinions, `*`))
# 
#   # Function to make frequency matrix
#   freq_matrix <- function(mat, levels_vec) {
#     sapply(seq_len(ncol(mat)), function(i) {
#       tab <- table(factor(mat[, i], levels = levels_vec))
#       as.numeric(tab)
#     })
#   }
# 
#   # Frequency matrices
#   #freq_mat0 <- freq_matrix(tb0, levels_vec)
#   freq_mat <- freq_matrix(tb, levels_vec)
#   #
#   #rownames(freq_mat0) <- levels_vec
#   #colnames(freq_mat0) <- colnames(tb0)
#   #rownames(freq_mat) <- levels_vec
#   #colnames(freq_mat) <- colnames(tb)
# 
#   # Plot
#   #par(mfrow = c(1, 2))  # adjust layout as needed
# 
#   #barplot(freq_mat0, names.arg = paste("group",1:ncol(bipartite0)),
#   #        col = my_palette,
#   #        main = "Initial Opinions", beside = FALSE)
# 
#   barplot(freq_mat,  names.arg = paste("group",1:ncol(bipartite)),
#           col = my_palette,
#           main = "Opinions", beside = FALSE)
# }

