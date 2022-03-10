#' Occupation time curve
#'
#' Computes the occupation time curve based on the raw activity data
#'
#' @param Xt a vector of raw activity data on a grid of time points
#' @param grid_widths a vector showing the number of time units (can be \eqn{< 1}) each grid point contains. This vector should have the same length as \code{Xt}.
#' @param as a vector showing the grid of activity levels over which the occupation time is evaluated.
#'
#' @return a vector containing the observed occupation time curve based on \code{Xt}
#'
#' @export
#'
#' @examples
mu_a_hat = function(Xt, grid_widths, as) {
  n_a = length(as)
  act_profile = sapply(1:n_a, FUN = function (a) {
    sum(grid_widths[Xt > as[a]])
  })
  return(act_profile)
}
