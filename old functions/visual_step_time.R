#' visual_step_time
#'
#' @docType package
#'
#' @author Victor Bernal \email{victor.arturo.bernal@gmail.com}
#'
#' @name visual_step_time
#'
#' @param opinion_history opinion_history
#' @param num_opinions num_opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_time
#' # visual_step_time(opinion_history, num_opinions)
#' #---------------------------------
#' @export
#visual_step_time <- function(opinion_history, num_opinions){
visual_step_time <- function(frac_mat, num_opinions){
  
  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  plot.new()
  #plot.window(xlim = c(1, ncol(opinion_history)), ylim = c(0, 1))
  plot.window(xlim = c(1, nrow(frac_mat)), ylim = c(0, 1))
  axis(1)
  axis(2)
  box()

  # frac <- apply(opinion_history, 2, function(x) prop.table(table(x)))
  # frac_mat <- matrix(unlist(frac), num_opinions, ncol(opinion_history) )

  # levels <- sort(unique(as.vector(opinion_history)))
  # n_ind <- nrow(opinion_history)
  # 
  # frac_mat <- sapply(levels, function(op) {
  #   colSums(opinion_history == op) / n_ind
  # })

  #frac_mat <- t(frac_mat)

  for (k in 1:num_opinions) {
    lines(x = 1:nrow(frac_mat), y = frac_mat[,k ],
          lwd = 2, col = my_palette[k])
  }

  # Update line plot for fraction of +1
  #Sys.sleep(0.05) # pause so you can see the update
}
