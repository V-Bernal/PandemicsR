#' Reconstruct a Bipartite Matrix from members
#'
#'
#' @author Victor Bernal
#'
#' @name reconstruct_bipartite
#'
#' @param members membership
#' @param n number of individuals
#' @param m number of groups
#'
#'
#' @return Bipartite Matrix.
#'
#' @examples
#' #---------------------------------
#' # reconstruct_bipartite(members, n, m)
#' #---------------------------------
#' @export
reconstruct_bipartite <- function(members, n, m) {
  i <- integer(0); j <- integer(0)
  for (g in seq_len(m)) {
    ids <- members[[g]]
    if (length(ids)>0) { i <- c(i,ids); j <- c(j, rep.int(g,length(ids))) }
  }
  sparseMatrix(i=i, j=j, x=1L, dims=c(n,m))
}
