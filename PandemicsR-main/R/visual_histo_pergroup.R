#' visual_histo_pergroup
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_histo_pergroup
#'
#' @import graphics
#'
#' @param opinions opinions
#' @param num_opinions number of opinions
#' @param members members
#'
#' @return visual_histo_pergroup \code{x}.
#'
#' @examples
#' #---------------------------------
#' # visual_histo_pergroup
#' # visual_histo_pergroup_multi(opinions, num_opinions, members)
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
