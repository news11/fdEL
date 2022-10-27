#' Cao and Cao2 confidence bands for functional means
#'
#' Compute the Cao and Cao2 confidence bands modified from the bands by Cao et al. (2012), as mentioned in Section 3, along with whether they capture the true functional mean or not.
#'
#' @param Y a matrix with the i-th row representing the i-th subject's observed functional data (e.g., occupation time data), the j-th column corresponding to the observed functional data at the j-th design time point (e.g., activity level) for each subject
#' @param as a vector consisting of the design time points of the observed functional data; that is, \eqn{{\bf G}_{n}} (from the smallest to the largest).
#' @param checkpts a vector consisting of the time points (from the smallest to the largest; e.g., the activity levels of interest) at which we want to evaluate the confidence band.
#' @param ub the largest time point (e.g., the activity levels of interest) at which we want to evaluate the confidence band.
#' @param mean_true a vector containing the true functional mean over \code{checkpts} (defined above)
#' @param c a parameter for controlling the number of interior knots in estimating the functional mean, as described in Section 4 of Cao et al. (2012) (default = 0.5)
#' @param p the order of spline functions (default = 2)
#' @param b_M the number of copies in simulating \eqn{\zeta(x)} for constructing the confidence bounds, as described in Section 4 of Cao et al. (2012) (default = 1000)
#' @param alpha a vector of significance levels of interest (default = c(0.05, 0.01))
#' @param N_m the number of interior knots in estimating the functional mean, as described in Section 4 of Cao et al. (2012)
#' @param cov_nnd 1 if the initially smoothed covariance estimates are projected onto the space of non-negative definite matrices and 0 if not; -1 if not doing the aforementioned projection and not taken into account the fact that the resulting estimated covariance can have complex eigenvalues
#'
#' @return a list containing the following:
#'   \item{cover_or_not}{a \code{length(alpha)}-vector with the \code{i}-th element indicating if the true functional mean is covered by the \eqn{100(1-}\code{alpha[i]}\eqn{)%} confidence band (1 if not covered and 0 otherwise)}
#'   \item{mean}{a vector containing the estimated mean function using spline estimation}
#'   \item{confidence_bounds}{a matrix containing the confidence bounds. The (\code{2 * i - 1})-th and (\code{2 * i})-th columns contain the upper and lower bound of the \eqn{100(1-}\code{alpha[i]}\eqn{)%} confidence band, respectively}
#'   \item{covariance}{the estimated covariance of the functional data; the initially smoothed covariance estimates are projected onto the space of non-negative definite matrices}
#'   \item{number_of_noises}{the number of eigenfunctions described in Section 4 of Cao et al. (2012)}
#'   \item{noises}{a matrix containing the \eqn{\hat{\phi}_{k}(\cdot)} functions defined in Section 4 of Cao et al. (2012), with the \eqn{k}-th column containing \eqn{\hat{\phi}_{k}(\cdot)}}
#'   \item{confidence_level}{\code{1 - alpha}}
#'   \item{runtime}{the runtime of the \code{Caoband} funciton}
#'
#' @importFrom splines bs
#' @importFrom stats lm
#' @importFrom stats rnorm
#' @importFrom stats quantile
#'
#' @export
#'
#' @references Cao, G., Yang, L. and Todem, D. (2012) Simultaneous inference for the mean function based on dense functional data. \emph{Journal of Nonparametric Statistics}, \strong{24}, 359--377.
#'
#' @author The code is first written by Yan-Yu Chen on 2021/02/03 and later modified slightly by Hsin-wen Chang.
#'
#' @examples
#' CaoCB <- Caoband(Y = OTdata$Ta[[1]], as = OTdata$as, checkpts = OTdata$as, mean_true = rep(0, 500))
#' CaoCB$confidence_bounds
Caoband <- function(Y, as, checkpts = rep(1:100)/100 * ub, ub = max(as), mean_true, c = 0.5, p = 2, b_M = 1000, alpha = c(0.05, 0.01), N_m = NA, cov_nnd = 1){
  start_time <- Sys.time()
  n <- dim(Y)[1]
  N <- dim(Y)[2]
  if (is.na(N_m) == 1) {
    N_m <- floor_n(c*n^(1 / (2 * p)) * log(n))
  }
  N_G <- floor_n(n ^ (1/(2 * p) * log(log(n))))
  X <- bs(as, knots = rep(1:N_m)/(N_m+1) * ub, degree = p-1, Boundary.knots = c(0, ub), intercept = TRUE)
  fit <- lm(colSums(Y)/n ~ X-1)
  beta_hat <- fit$coefficients
  m_hat_check <- bs(checkpts, knots = rep(1:N_m) / (N_m+1) * ub, degree = p - 1, Boundary.knots = c(0, ub), intercept = TRUE) %*% beta_hat
  m_hat_matrix <- matrix(rep(1,n * N), ncol = N, nrow = n) %*% diag(fit$fitted.values)
  C_matrix <- t(Y - m_hat_matrix) %*% (Y - m_hat_matrix)/n
  C <- c(C_matrix[col(C_matrix) != row(C_matrix)])
  G_matrix <- bs(as, knots = rep(1:N_G) / (N_G + 1) * ub, degree = p - 1, Boundary.knots = c(0, ub), intercept = TRUE)
  B <- matrix(rep(0,(N ^ 2 - N) * (N_G + p) ^ 2), nrow = N ^ 2 - N, ncol = (N_G + p) ^ 2)
  idx <- 1
  for (i in 1:(N_G + p)){
    for (j in 1:(N_G + p)){
      V <- G_matrix[, j] %*% t(G_matrix[, i])
      B[, idx] = c(V[col(V) != row(V)])
      idx = idx + 1
    }
  }
  b <- matrix(lm(C ~ B - 1)$coefficients, ncol = N_G + p, nrow = N_G + p)
  checkpts_matrix <- bs(checkpts, knots = rep(1:N_G)/(N_G + 1) * ub, degree = p - 1, Boundary.knots = c(0, ub), intercept = TRUE)
  if (cov_nnd == 1) {
    G_hat_0 <- checkpts_matrix %*% b %*% t(checkpts_matrix)
    eigen_G_0 <- eigen(t(G_hat_0))
    T <- which(Re(eigen_G_0$values) <= 0)[1] - 1
    if (T == 1) {
      G_hat = as.matrix(eigen_G_0$vectors[, 1:T]) %*% as.matrix(diag(as.matrix(eigen_G_0$values[1:T]))) %*% t(as.matrix(eigen_G_0$vectors[, 1:T]))
    } else {
      G_hat = as.matrix(eigen_G_0$vectors[,1:T]) %*% as.matrix(diag(eigen_G_0$values[1:T])) %*% t(as.matrix(eigen_G_0$vectors[,1:T]))
    }
    G_hat = Re(G_hat)
    eigen_G <- eigen(t(G_hat))
  } else if (cov_nnd %in% c(0, -1)) {
    G_hat <- checkpts_matrix %*% b %*% t(checkpts_matrix)
    if (cov_nnd == 0) {
      G_hat = Re(G_hat)
    }
    eigen_G <- eigen(t(G_hat))
    if (cov_nnd == 0) {
      T <- which(Re(eigen_G$values) <= 0)[1] - 1
    } else if (cov_nnd == -1) {
      T <- which(eigen_G$values <= 0)[1] - 1
    }
  }
  if (cov_nnd %in% c(0, 1)) {
    kappa <- which(cumsum(Re(eigen_G$values)[1:T])/sum(Re(eigen_G$values)[1:T]) > 0.95)[1] # number of noises
    if (kappa == 1){
      phi_hat <- Re(eigen_G$vectors)[,1:kappa] * sqrt(Re(eigen_G$values)[1:kappa])
    } else{
      phi_hat <- Re(eigen_G$vectors)[, 1:kappa] %*% diag(c(sqrt(Re(eigen_G$values)[1:kappa])))
    }
    phi_hat = Re(phi_hat)
  } else if (cov_nnd == -1) {
    kappa <- which(cumsum(eigen_G$values[1:T])/sum(eigen_G$values[1:T]) > 0.95)[1] # number of noises
    if (kappa == 1){
      phi_hat <- eigen_G$vectors[,1:kappa] * sqrt(eigen_G$values[1:kappa])
    } else{
      phi_hat <- eigen_G$vectors[, 1:kappa] %*% diag(c(sqrt(eigen_G$values[1:kappa])))
    }
  }
  xi_hat_m <- apply(abs(sqrt(diag(1 / diag(G_hat))) %*% phi_hat %*% matrix(rnorm(b_M * kappa, 0, 1), ncol = b_M, nrow = kappa)), 2, max)
  Q_hat <- quantile(xi_hat_m, probs = 1 - alpha, names = FALSE)
  range_CB <- sqrt(diag(G_hat)) %*% t(Q_hat) / sqrt(n)
  CB <- matrix(rep(0, 2 * length(alpha) * length(checkpts)), nrow = length(checkpts), ncol = 2 * length(alpha))
  cover_or_not <- rep(0, length(alpha))
  col_names <- rep(0, 2 * length(alpha))
  for (i in rep(1:length(alpha))){
    CB[, 2 * i - 1] = m_hat_check + range_CB[, i]
    CB[, 2 * i] = m_hat_check - range_CB[, i]
    col_names[2 * i - 1] = paste('upper_', toString(1 - alpha[i]), sep = '')
    col_names[2 * i] = paste('lower_', toString(1 - alpha[i]), sep = '')
    if (any(CB[, 2 * i - 1] < mean_true) || any(CB[, 2 * i] > mean_true)){
      cover_or_not[i] = 1
    }
  }
  colnames(CB) <- col_names
  end_time <- Sys.time()
  res <- list('cover_or_not' = cover_or_not, 'mean' = m_hat_check, 'confidence_bounds' = CB, 'covariance' = G_hat, 'number_of_noises' = kappa, 'noise' = phi_hat, 'confidence_level' = 1-alpha, 'runtime' = end_time-start_time)
  return(res)
}
