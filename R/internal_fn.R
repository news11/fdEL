#' @importFrom stats median
#' @importFrom stats var
#' @importFrom nloptr nloptr
#' @importFrom methods is

floor_n = function(x, n = 1) {
  # Extending the floor function so that for the given input number x and positive integer n,
  # the function returns the number b / n for some interger b, such that b / n <= x < (b + 1) / n
  # The round function is to eliminate a small numerical difference from the actual value of x * n.
  # Args:
  #   x : a real number
  #   n: a positive integer n
  # Returns:
  #   the number b / n for some interger b, such that b / n <= x < (b + 1) / n
  floor(round(x * n, 2)) / n
}
product_mat <- function(x, y) {
  # A modified * function so that  0 * Inf = 0
  # Args:
  #   x, y: two numbers or vectors to be multiplied
  # Returns:
  #   x * y, where an element becomes 0 if the corresponding element of one of x or y is 0
  out <- x*y
  for (i in 1:dim(x)[1]) {
      out[i, (x[i, ] == 0 | y[i, ] == 0)] <- 0
  }
  return (out)
}
product_arr3 <- function(x, y) {
  # A modified product_mat function so that  0 * Inf = 0, the 3 dimensional array version
  # Args:
  #   x, y: two 3 dimensional arrays to be multiplied
  # Returns:
  #   x * y, where an element becomes 0 if the corresponding element of one of x or y is 0
  out <- array(0, dim(x))
  for (j in 1:dim(x)[3]) {
    out[, , j] = product_mat(x[, , j], y[, , j])
  }
  return (out)
}
U_division00_mat = function(x, y) {
  # A modified product_mat function so that  0 * Inf = 0, the matrix version
  # Args:
  #   x, y are matrices of the same dimensions
  # Returns:
  #   out: x * y; when x[i, ] = 0 or y[i, ] = 0, out[i, ] = 0
  out = x / y
  for (i in 1:dim(x)[1]) {
    out[i, (x[i, ] == 0 & y[i, ] == 0)] = 0
  }
  return(out)
}
log_likelihood_p1_lamb <- function(lbs, m_lbs, init_lamb, maxeval = 10000) {
  # Computes -1 times the numerator of the local EL ratio at a fixed design time point \eqn{a \in {\bf G}_n}, for constructing the EL confidence band
  # Args:
  #   lbs: an n-vector \{f_n(T_{1})(a), .., f_n(T_{n_j,j})(a))
  #   m_lbs: a number: \tilde{\mu}(a)
  #   init_lamb: an initival value for \tilde{\lambda}(a) in optimization of R(\tilde{\mu})(a); see Supplement Section 4 for the definition of \tilde{\lambda}(a)
  #   maxeval: the termination condition for evaluating the objective function (i.e., -1 times the log likelihood) in the optimization of the numerator of R(\tilde{\mu})(a):
  #   the maximum number of objective function evaluations
  # Returns:
  #   nlopt: a list (some of whose elements are described in the following lines; the rest of the elements are from the return value of the function nloptr, see \url{https://CRAN.R-project.org/package=nloptr})
  #   of outputs from optimizing the numerator of R(\tilde{\mu})(a) using the function nloptr from the R package nloptr:
  #   solution: the resulting \tilde{\lambda}(a) defined in Supplement Section 4
  #   objective: the objective function -\sum_{i=1}^n \log \tilde{p}_i(a), where \tilde{p}_i(a) is defined in Supplement Section 4
  #   status: an integer value with the status of the optimization, fully described in \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values};
  #   we care about the value 5 when the optimization is stopped only because maxeval (above) was reached, but the optimization is deemed successful and returns no error in R
  #   resultEEs: the value of the constraint \sum_{i=1}^n \tilde{p}_i(a) \tilde{Y}_i(a), where \tilde{Y}_i(a) is defined in Supplement Section 4
  #   p: the n-vector (\tilde{p}_1(a), ..., \tilde{p}_n(a))
  h <- length(lbs)
  fnc_objective <- function(lamb) {
    p = 1 / (h * (1 + lamb * (lbs - m_lbs)))
    x <- log(p)
    return (-sum(x))
  }
  gnc_objective <- function(lamb) {
    p = 1 / (h * (1 + lamb * (lbs - m_lbs)))
    return (sum(h * p * (lbs - m_lbs)))
  }
  constraints <- function(lamb) {
    p = 1 / (h * (1 + lamb * (lbs - m_lbs)))
    return(sum(p * (lbs - m_lbs)))
  }
  jac_constraints <- function(lamb) {
    p = 1 / (h * (1 + lamb * (lbs - m_lbs)))
    return(-sum(h * (p ^ 2) * (lbs - m_lbs) ^ 2))
  }
  local_opts <- list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-10)
  opts_used  <- list("algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-10, "maxeval" = maxeval, "local_opts" = local_opts)
  ll_l = max((1 / h - 1) / (lbs - m_lbs)[(lbs - m_lbs) > 0])
  ul_l = min((1 / h - 1) / (lbs - m_lbs)[(lbs - m_lbs) < 0])
  nlopt <- nloptr(x0 = init_lamb, eval_f = fnc_objective, eval_grad_f = gnc_objective,
                  lb = ll_l, ub = ul_l,
                  eval_g_eq = constraints, eval_jac_g_eq = jac_constraints, opts = opts_used)
  nlopt$resultEEs = constraints(nlopt$solution)
  nlopt$p = 1 / (h * (1 + nlopt$solution * (lbs - m_lbs)))
  return(nlopt)
}
model = function(lamb_gam, lbs, n_subject) {
  # Computes the left hand side of the estimating equations needed to obtain
  # pointwise EL statistic for k samples at a fixed design time point \eqn{a \in {\bf G}_n}
  # Args:
  #   lamb_gam: a vector (lambda_j(a), j= 1, ..., k-1, gamma_j(a), 1, ..., k), where gamma_j(a) corresponds to \gamma_j + {lambda_{j-1}(a) - lambda_{j}(a)} * {-\hat{\mu}_j(a)} defined in Supplement Section 8
  #   lbs: for fixed a, a list whose j-th entry is the n_j-vector (T_{1j}(a), .., T_{n_j,j}(a)) in Supplement Section 8
  #   n_subject: a k-vector of sample size from each group x
  # Returns:
  #   out: (\sum_i p_{i, j + 1}(a) T_{i, j + 1}(a) - \sum_i p_{ij}(a) T_{i, j}(a), j = 1, .., n_group - 1,
  #        \sum_i p_{ij}(a) - 1, j = 1, ..., n_group)
  n = sum(n_subject)
  gamma_js = n_subject / n
  n_group = length(n_subject)
  out = rep(0, times = 2 * n_group - 1)
  lamb = lamb_gam[1:(n_group - 1)]
  gam = lamb_gam[-(1:(n_group - 1))]
  lamb_kp1 = c(0, lamb, 0)
  p = list()
  p[[1]] = 1 / (n * (gam[1] + (lamb_kp1[1] - lamb_kp1[2]) * lbs[[1]]))
  EEs = sapply(1:n_group, FUN = function (j) {
    out[n_group + j - 1] <<- sum(p[[j]]) - 1
    if (j != n_group){
      p[[j + 1]] <<- 1 / (n * (gam[j + 1] + (lamb_kp1[j + 1] - lamb_kp1[j + 2]) * lbs[[j + 1]]))
      out[j] <<- sum(p[[j + 1]] * lbs[[j + 1]]) - sum(p[[j]] * lbs[[j]])
    }
  })
  return(out)
}
neg2logR_CB = function(mu_tilde, Tai, EL_CBcrit = 0, error_out = 0, parameters = list(n_lamb_grid = 1000, tol_resultEEs = 10 ^ (-4))) {
  # Computes -2log R(\tilde{\mu})(a) for a fixed design time point \eqn{a \in {\bf G}_n}
  # Args:
  #   mu_tilde: a vector \tilde{\mu}(a) of length length(as), to be solved for in confidence band construction
  #   Tai: a vector, containing the observed functional data (e.g., occupation time data) at the design time point (e.g., activity level) \eqn{a} for each subject
  #   EL_CBcrit: critical value for constructing the confidence band
  #   parameters: a list (whose elements are described in the following lines) of parameters related to initilization and tolerance in optimization of R(\tilde{\mu})(a):
  #   n_lamb_grid: number of grid points for initializing lambda(a), the Lagrange multiplier related to optimization of the numerator of R(\tilde{\mu})(a); see Supplement Section 4)
  #   tol_resultEEs: tolerance for the value of the contraints in the optimization of the numerator of R(\tilde{\mu})(a)
  # Returns:
  #   -2log R(\tilde{\mu})(a)
  n_subject = length(Tai)
  ll_l = max((1 / n_subject - 1) / (Tai - mu_tilde)[(Tai - mu_tilde) > 0])
  ul_l = min((1 / n_subject - 1) / (Tai - mu_tilde)[(Tai - mu_tilde) < 0])
  init_lambda = median(c(ll_l, ul_l))
  num_neg_loglik = log_likelihood_p1_lamb(lbs = Tai, m_lbs = mu_tilde, init_lamb = init_lambda)
  denom_neg_loglik = n_subject * log(n_subject)
  num_warn <- tryCatch(log_likelihood_p1_lamb(lbs = Tai, m_lbs = mu_tilde, init_lamb = init_lambda), error = function(e) e, warning = function(w) w)
  init_lambda_grid = seq(ll_l, ul_l, length = parameters$n_lamb_grid)
  init_ind = 1
  status_grid = matrix(0, nrow = parameters$n_lamb_grid, ncol = 5)
  while ((abs(num_neg_loglik$resultEEs) >= parameters$tol_resultEEs) != 0 | is(num_warn, "warning") != 0 | num_neg_loglik$status == 5) {
    if (init_ind > parameters$n_lamb_grid) break  # break out of the while loop
    init_lambda = init_lambda_grid[init_ind]
    num_neg_loglik = log_likelihood_p1_lamb(lbs = Tai, m_lbs = mu_tilde, init_lamb = init_lambda)
    num_warn <- tryCatch(log_likelihood_p1_lamb(lbs = Tai, m_lbs = mu_tilde, init_lamb = init_lambda), error = function(e) e, warning = function(w) w)
    status_grid[init_ind, ] = c(num_neg_loglik$status, num_neg_loglik$resultEEs, num_neg_loglik$solution, num_neg_loglik$objective, sum(num_neg_loglik$p))
    init_ind = init_ind + 1
  }
  still_error = ((abs(num_neg_loglik$resultEEs) >= parameters$tol_resultEEs) != 0 | is(num_warn,"warning") != 0 | num_neg_loglik$status == 5)
  test1t = 2 * num_neg_loglik$objective - 2 * denom_neg_loglik
  if (error_out == 0) {
    return(test1t - EL_CBcrit)
  } else {
    return(c(test1t - EL_CBcrit, still_error))
  }
}
neg2logR_test = function(ai, Ta, n_subject, EL_testcrit = 0, error_out = 0, parameters = list(n_mu_grid = 20, n_lamb_grid = 20, n_grid_incbd = 6, mu_tol_grid = 20, lamb_tol_grid = 20)) {
  # Computes -2logR_k(a) - EL_testcrit for a fixed design time point \eqn{a \in {\bf G}_n}
  # Args:
  #   ai: index for the given design time point \eqn{a}
  #   Ta: a list with the j-th element being the matrix of observed functional data (e.g., occuptation time curve), i-th row representing the i-th subject's observed curve
  #   (dim = n_subject[j] x length(as[[j]]))
  #   n_subject: a vector, n_subject[j] being the sample size for the j-th group
  #   EL_testcrit: critical value for testing H0
  #   error_out: = 1 to output more information about error report (if any); = 0 to output -2logR_k(a) - EL_testcrit only
  #   parameters: a list (whose elements are described in the following lines) of parameters related to initilization and tolerance in optimization of R_k(a):
  #   n_mu_grid: number of grid points for mu initialization,
  #   n_lamb_grid: number of grid points for lambdas initialization
  #   n_grid_incbd: max number of times the number of grid points can be updated
  #   mu_tol_grid: a parameter that makes mu initialization away from max of mininums in each group, and min of maximums in each group
  #   lamb_tol_grid: a parameter that makes lambdas initialization away from lower and upper bounds of lambdas
  # Returns:
  #   -2logR_k(a) - EL_testcrit
  # Returns only if error_out = 1:
  #   still_error: = 1 if there is error in optimization in R_k(a) computing; 0 otherwise
  #   status = c(whether optimization in R_k(a) returns meaningful values, precision in optimization,
  #            the value of the output from the function `model' evaluated at `root' below, a list of k-elements with j-th element being minimal p_{ij}(a) for group j)
  #   roots = root from solving estimating equation for optimization in R_k(a)
  n = sum(n_subject)
  gamma_js = n_subject / sum(n_subject)
  n_group = length(n_subject)
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range((Ta[[j]])[, ai])
  }), nrow = 2)
  # dealing with non-overlapping regions (include the case when the region is degenerate, i.e., the = case)
  if (max(min_ups_max_lows[1, ]) >= min(min_ups_max_lows[2, ])) {
    if (error_out == 0) {
      return(Inf)
    } else {
      return(list(neg2logR_crit = Inf,
        still_error = 0
      ))
    }
  }
  # END dealing with non-overlapping regions
  T_js = lapply(1:n_group, FUN = function (j) {
     (Ta[[j]])[, ai] / (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ]))
  })
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range(T_js[[j]])
  }), nrow = 2)
  n_grid_inc = 0
  still_error = 1
  test1t = Inf
  n_mu_grid = parameters$n_mu_grid
  n_lamb_grid = parameters$n_lamb_grid
  while (still_error == 1 & n_grid_inc <= parameters$n_grid_incbd) {
    n_grid_inc = n_grid_inc + 1
    n_mu_grid = n_mu_grid * n_grid_inc + 1
    n_lamb_grid = n_lamb_grid * n_grid_inc + 1
    # grid construction
    mu_tol = (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ])) / parameters$mu_tol_grid
    init_mu_a_grid = seq(max(min_ups_max_lows[1, ]) + mu_tol, min(min_ups_max_lows[2, ]) - mu_tol, length = n_mu_grid)
    ll_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    ul_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    init_Delta_grid = array(unlist(sapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      sapply(1:n_group, FUN = function (j) {
        g_j = (T_js[[j]] - init_mu_a) / gamma_js[j]  # a vector of length n_j
        ll_l[init_mu_ind, j] <<- max((1/n_subject[j] - 1) / g_j[g_j > 0])
        ul_l[init_mu_ind, j] <<- min((1/n_subject[j] - 1) / g_j[g_j < 0])
        lamb_tol = (ul_l[init_mu_ind, j] - ll_l[init_mu_ind, j]) / parameters$lamb_tol_grid
        seq(ll_l[init_mu_ind, j] + lamb_tol, ul_l[init_mu_ind, j] - lamb_tol, length = n_lamb_grid)
      })
    })), c(n_lamb_grid, n_group, n_mu_grid))
    # END grid construction

    # initialization
    init_check = matrix(0, nrow = n_mu_grid, ncol = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
    init_Delta_grid_expand = matrix(unlist(lapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      mat1 = cbind(init_mu_a, sapply(1:(n_group - 1), FUN = function (j) {
        rep(rep(init_Delta_grid[, j, init_mu_ind],each = n_lamb_grid ^ (j - 1)), length = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
      }))

      init_lamb_km1 = -apply(as.matrix(mat1[, -1]), 1, sum)
      init_check[init_mu_ind, ] <<- (init_lamb_km1 > ll_l[init_mu_ind, n_group] & init_lamb_km1 < ul_l[init_mu_ind, n_group])
      return(t(mat1[(init_check[init_mu_ind, ] ==1), ]))
    })), byrow = T, ncol = n_group)
    # END initialization
    n_init = dim(init_Delta_grid_expand)[1]
    status_grid = matrix(NA, nrow = n_init, ncol = 2 + 3 * n_group - 1)
    roots = matrix(NA, nrow = n_init, ncol = 2 * n_group - 1)
    multir = list()
    multir$estim.precis = 1
    min_ps = rep(-1, times = n_group)
    for (init_ind in 1:n_init) {
      init_mu_a = init_Delta_grid_expand[init_ind, 1]
      init_lambda_js = cumsum(-init_Delta_grid_expand[init_ind, -1])
      init_gam = gamma_js - diff(c(0, init_lambda_js, 0)) * (-init_mu_a)
      lamb_gam = c(init_lambda_js, init_gam)
      options(warn = 2)
      possibleError <- tryCatch(
        multir<-rootSolve::multiroot(f = model, start = lamb_gam, lbs = T_js, n_subject = n_subject), #can input multir here!
        error=function(e) e
      )
      options(warn = 0)
      if (!inherits(possibleError, "error")) {
        roots[init_ind, ] = multir$root
        lamb = roots[init_ind, ][1:(n_group - 1)]
        gam = roots[init_ind, ][-(1:(n_group - 1))]
        lamb_kp1 = c(0, lamb, 0)
        p = list()
        min_ps = sapply(1:n_group, FUN = function (j) {
          p[[j]] <<- 1 / (n * (gam[j] + (lamb_kp1[j] - lamb_kp1[j + 1]) * T_js[[j]]))
          min(p[[j]])
        })
        status_grid[init_ind, ] = c(!inherits(possibleError, "error"), multir$estim.precis, multir$f.root, min_ps)
        if (multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0) break
      }
    }
    still_error = !(multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0)
  }
  if (exists("p") == 1) {
    test1t = -2 * sum(unlist(lapply(1:n_group, FUN = function (j) {
      return(log(p[[j]]))
    }))) - 2 * sum(n_subject * log(n_subject))
  }
  if (error_out == 0) {
    return(test1t - EL_testcrit)
  } else {
    return(list(neg2logR_crit = test1t - EL_testcrit,
      still_error = still_error,
      status = status_grid[init_ind, ],
      roots = roots[init_ind, ])
    )
  }
}
