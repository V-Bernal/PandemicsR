#' visual_step_time
#'
#' @author Victor Bernal
#'
#' @name visual_step_time
#'
#' @import graphics
#'
#' @param frac_mat fraction of opinions
#' @param num_opinions number of opinions
#'
#' @return NULL
#'
#' @examples
#' #---------------------------------
#' # visual_step_time(frac_mat, num_opinions)
#' #---------------------------------
#' @export
visual_step_time <- function(frac_mat, num_opinions){

  # increase top margin
  op <- par(mar = c(5, 4, 6, 2))  # bottom, left, top, right


  levels_vec <- get_levels_vec(num_opinions)
  my_palette <- get_palette(num_opinions)

  # extract time column
  time_vals <- frac_mat[,1]

  # opinion fraction columns
  frac_vals <- frac_mat[, -1, drop = FALSE]

  plot.new()

  plot.window(
    xlim =  c(min(time_vals), max(time_vals)) ,
    ylim = c(0,1),
    xaxt = "n", bty='L')

  # evenly spaced tick marks
  ticks <- pretty(time_vals, n = 6)

  axis(1,
       at = ticks,
       labels = round(ticks, 2))
  axis(2)

  box()

  title(
    xlab = "Time",
    ylab = "Voter Fraction"
  )

  for (k in seq_len(ncol(frac_vals))) {

    lines(
      x = time_vals,
      y = frac_vals[,k],
      type = "s",
      lwd = 2,
      col = my_palette[k]
    )
  }
}
