#' EL and related confidence bands for functional means
#'
#' Computes EL-based confidence bands and related Wald-type confidence bands for the mean of the functional data of interest, along with whether they capture the true functional mean or not
#'
#' @param Ta a matrix with the i-th row representing the i-th subject's observed functional data (e.g., occupation time data), the j-th column corresponding to the observed functional data at the j-th design time point (e.g., activity level) for each subject
#' @param as a vector consisting of the design time points of the observed functional data; that is, \eqn{{\bf G}_n} defined in Section 2.2 (from the smallest to the largest).
#' @param as_eval a vector consisting of the time points (from the smallest to the largest; e.g., the activity levels of interest) at which we want to evaluate the confidence band. This vector can be denser or coarser than \code{as}, as long as the largest point of \code{as_eval} is no greater than the last point of \code{as} (defined above).
#' @param n_boot number of bootstrap samples.
#' @param mu_a a vector containing the true functional mean over \code{as_eval} (defined above)
#' @param alpha a number, the significance level of interest
#' @param as_right_ind_keep0 This is the index in the vector \code{as_eval} (defined above) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1. The default value is \code{length(as_eval)}, where the approach in Supplement Section 5.1 is NOT implemented.
#' @param as_right_ind_keep_as0 This is the index in the vector \code{as} (defined above) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1. The default value is \code{length(as)}, where the approach in Supplement Section 5.1 is NOT implemented.
#' @return a list containing the following:
#'   \item{EL_CB_Nair}{EL confidence band over \code{as_eval} (defined above) after adapting Nair's two-step approach described in Supplement Section 5.1, with first row being the lower bound and second row being the upper bound}
#'   \item{EL_CB_Nairrej}{1 if the truth \code{mu_a} is not included in the EL confidence band; 0 otherwise}
#'   \item{EL_CBcrit}{critical value for computing the EL confidence band}
#'   \item{EP_CB_Nair}{EP confidence band over \code{as_eval} (defined above) after adapting Nair's two-step approach described in Supplement Section 5.1, with first row being the lower bound and second row being the upper bound}
#'   \item{EP_CB_Nairrej}{1 if the truth \code{mu_a} is not included in the EP confidence band; 0 otherwise}
#'   \item{EP_CBcrit}{critical value for computing the EP confidence band}
#'   \item{HW_CB_eval}{NS confidence band over \code{as_eval} (defined above), with first row being the lower bound and second row being the upper bound}
#'   \item{HW_CBrej}{1 if the truth \code{mu_a} is not included in the NS confidence band; 0 otherwise}
#'   \item{HW_CBcrit}{critical value for computing the NS confidence band}
#'   \item{mu_hat_vec}{a \code{length(as)}-vector containing \eqn{\hat{\mu}(a)}, the sample mean of the observed functional data}
#'   \item{S2_hat_vec}{a \code{length(as)}-vector containing \eqn{\hat{S}^2(a)}, the estimated variance of the observed functional data defined in Section 2.3}
#'   \item{EL_CB_aicheck}{an array, whose (bi, ai, zi) entry recording the following information for the bi-th bound (bi = 1 for lower bound, bi = 2 for upper bound) of the EL band at \code{as[ai]} when adapting Nair's two-step approach described in Supplement Section 5.1 using the zi-th element in the vector of \code{as_right_ind_keep_as0} (defined above): the entry is 1 if the bound EL band is no greater than the minimum of the observed functional data at \code{as[ai]} plus some tolerance value, or when the bound of the EL band is too close to the extreme values of the observed functional data at \code{as[ai]}; 0 otherwise}
#'   \item{neg2logR_CB_decisionstop_error}{1 if there is an error during the EL confidence band construction; 0 otherwise}
#'
#' @importFrom stats sd
#' @importFrom stats uniroot
#'
#' @export
#'
#' @examples
#' CB_out <- elfband(Ta = OTdata$Ta[[1]], as = OTdata$as, as_eval = OTdata$as)
#' CB_out$EL_CB_Nair
elfband = function(Ta, as, as_eval, n_boot = 1000, mu_a = NA, alpha = 0.05, as_right_ind_keep0 = length(as_eval), as_right_ind_keep_as0 = length(as)) {
  n_Ta = length(as)
  n_subject = dim(Ta)[1]
  n_agrid = length(as)
  n_agrid_eval = length(as_eval)
  right_endpt_len = length(as_right_ind_keep_as0)
  boot_indx = sample(1:n_subject, n_subject * n_boot, replace = T)
  boot_indx_mat = matrix(boot_indx, byrow = T, nrow = n_boot, ncol = n_subject)
  mu_hat = matrix(rep(apply(Ta, 2, mean), each = n_boot), nrow = n_boot, ncol = n_Ta)
  Ta_boot = array(Ta[boot_indx, ], c(n_subject, n_boot, n_Ta))
  mu_hat_star = apply(Ta_boot, c(2, 3) , mean)
  S2_hat = matrix(apply(Ta, 2, var) * (n_subject - 1) / n_subject, byrow = T, nrow = n_boot, ncol = n_Ta)
  U_hat_star = sqrt(n_subject) * U_division00_mat(mu_hat_star -  mu_hat, sqrt(S2_hat))
  U2_hat_star = U_hat_star ^ 2
  b_a_vec = t(sapply(1:n_agrid_eval, FUN = function (i) {
    min(which(as >= as_eval[i]))
  }))
  EL_CBdistr = list()
  EL_CBcrit = 1:right_endpt_len * 0
  EP_CBdistr = list()
  EP_CBcrit = 1:right_endpt_len * 0
  EP_CB = array(-1, c(2, n_Ta, right_endpt_len))
  EP_CB_eval = array(-1, c(2, n_agrid_eval, right_endpt_len))
  HW_CBdistr = list()
  HW_CBcrit = 1:right_endpt_len * 0
  HW_CB = array(-1, c(2, n_Ta, right_endpt_len))
  HW_CB_eval = array(-1, c(2, n_agrid_eval, right_endpt_len))
  for (as_right_ind in 1:right_endpt_len) {
    EL_CBdistr[[as_right_ind]] = apply(U2_hat_star[, 1:as_right_ind_keep_as0[as_right_ind]], 1, max)
    EL_CBcrit[as_right_ind] = as.vector(quantile(EL_CBdistr[[as_right_ind]], probs = 1 - alpha[1]))
    EP_CBdistr[[as_right_ind]] = apply(abs(U_hat_star[, 1:as_right_ind_keep_as0[as_right_ind]]), 1, max)
    EP_CBcrit[as_right_ind] = as.vector(quantile(EP_CBdistr[[as_right_ind]], probs = 1 - alpha[1]))
    EP_CB[1, , as_right_ind] = mu_hat[1, ] - EP_CBcrit[as_right_ind] / sqrt(n_subject) * sqrt(S2_hat[1, ])
    EP_CB[2, , as_right_ind] = mu_hat[1, ] + EP_CBcrit[as_right_ind] / sqrt(n_subject) * sqrt(S2_hat[1, ])
    EP_CB_eval[1, , as_right_ind] = (EP_CB[1, , as_right_ind])[b_a_vec]
    EP_CB_eval[2, , as_right_ind] = (EP_CB[2, , as_right_ind])[b_a_vec]
    HW_CBdistr[[as_right_ind]] = apply(abs(U_hat_star * sqrt(S2_hat))[, 1:as_right_ind_keep_as0[as_right_ind]], 1, max)
    HW_CBcrit[as_right_ind] = as.vector(quantile(HW_CBdistr[[as_right_ind]], probs = 1 - alpha[1]))
    HW_CB[1, , as_right_ind] = mu_hat[1, ] - HW_CBcrit[as_right_ind] / sqrt(n_subject)
    HW_CB[2, , as_right_ind] = mu_hat[1, ] + HW_CBcrit[as_right_ind] / sqrt(n_subject)
    HW_CB_eval[1, , as_right_ind] = (HW_CB[1, , as_right_ind])[b_a_vec]
    HW_CB_eval[2, , as_right_ind] = (HW_CB[2, , as_right_ind])[b_a_vec]
  }
  EL_CB = array(-1, c(2, n_Ta, right_endpt_len))
  EL_CB_aicheck = array(-1, c(2, n_Ta, right_endpt_len))
  tolfac_ls = matrix(0, nrow = n_Ta, ncol = right_endpt_len)
  tolfac_us = matrix(0, nrow = n_Ta, ncol = right_endpt_len)
  neg2logR_CB_decisionstop_error = 0
  prev2 = NA
  prev1 = NA
  m1m2lbub_digit = 6
  for (ai in 1:max(as_right_ind_keep_as0)) {
    a = as[ai]
    if (var(Ta[,ai]) == 0) {
      EL_CB[, ai, ] = rep(Ta[1, ai], 2 * right_endpt_len)
      prev2 = NA
      prev1 = NA
      next
    }
    if (ai > 1 & sum(Ta[, ai] != Ta[, (ai - 1)]) == 0) {
      EL_CB[, ai, ] = EL_CB[, (ai - 1), ]
      tolfac_ls[ai, ] = tolfac_ls[ai - 1, ]
      tolfac_us[ai, ] = tolfac_us[ai - 1, ]
      next
    }
    cen = mean(Ta[, ai])
    sc = sd(Ta[, ai])
    Tai = (Ta[, ai] - cen) / sc
    mu_a_hat = (mu_hat[1, ai] - cen) / sc
    tolfac_ai = 4
    tolfac_l = 10
    tolfac_u = 10
    min_a = min(Tai)
    mu_a_hat_orig = mu_a_hat
    max_a_orig = max(Tai)
    lb = min_a + 10 ^ (-tolfac_ai) * (mu_a_hat_orig - min_a)
    for (as_right_ind in 1:right_endpt_len) {
      if (ai == 1 | is.na(prev1)) {
        m1_start = mu_a_hat_orig
      } else if (((EL_CB[1, ai - 1, as_right_ind] - cen) / sc) > min_a & ((EL_CB[1, ai - 1, as_right_ind] - cen) / sc) > lb){
        mu_a_hat_alt = (EL_CB[1, ai - 1, as_right_ind] - cen) / sc + 2 * 10 ^ (-tolfac_ai) * (mu_a_hat_orig - min_a)
        m1_start = min(mu_a_hat_orig, mu_a_hat_alt)
      } else {
        m1_start = mu_a_hat_orig
      }
      m1 = round(m1_start - 10 ^ (-tolfac_ai) * (m1_start - min_a), m1m2lbub_digit)
      m2 = round(mu_a_hat + 10 ^ (-tolfac_ai) * (max_a_orig - mu_a_hat), m1m2lbub_digit)
      if (ai == 1 | is.na(prev2)) {
        max_a = max_a_orig
      } else if (((EL_CB[2, ai - 1, as_right_ind] - cen) / sc) < max_a_orig - 10 ^ (-tolfac_ai) * (max_a_orig - mu_a_hat) & ((EL_CB[2, ai - 1, as_right_ind] - cen) / sc) > m2) {
        max_a_alt = (EL_CB[2, ai - 1, as_right_ind] - cen) / sc + 10 ^ (-tolfac_ai) * (max_a_orig - mu_a_hat)
        max_a = min(max_a_orig, max_a_alt)
      } else {
        max_a = max_a_orig
      }
      ub = max_a - 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
      while (sum((Tai - m1) > 0) %in% c(0, n_subject) | sum((Tai - m2) > 0) %in% c(0, n_subject) | sum((Tai - lb) > 0) %in% c(0, n_subject) | sum((Tai - ub) > 0) %in% c(0, n_subject)) {
        tolfac_ai = tolfac_ai + 1
        m1 = round(m1_start - 10 ^ (-tolfac_ai) * (m1_start - min_a), m1m2lbub_digit)
        m2 = round(mu_a_hat + 10 ^ (-tolfac_ai) * (max_a - mu_a_hat), m1m2lbub_digit)
        ub = max_a - 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
        lb = min_a + 10 ^ (-tolfac_ai) * (m1_start - min_a)
      }
      possibleError11 <- tryCatch(
        testm1 <- neg2logR_CB(m1, Tai, EL_CBcrit[as_right_ind], error_out = 1),
        error = function(e) e
      )
      possibleError12 <- tryCatch(
        testm2 <- neg2logR_CB(m2, Tai, EL_CBcrit[as_right_ind], error_out = 1),
        error = function(e) e
      )
      possibleError2 <- tryCatch(
        testub <- neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind], error_out = 1),
        error = function(e) e
      )
      possibleError3 <- tryCatch(
        testlb <- neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind], error_out = 1),
        error = function(e) e
      )
      if((testm1[2] + testm2[2] + testub[2] + testlb[2] == 0) & !inherits(possibleError11, "error") & !inherits(possibleError12, "error") & !inherits(possibleError2, "error") & !inherits(possibleError3, "error")) {
        testm1 = neg2logR_CB(m1, Tai, EL_CBcrit[as_right_ind])
        testm2 = neg2logR_CB(m2, Tai, EL_CBcrit[as_right_ind])
        testub = neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind])
        testlb = neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind])
        tolfac_ai1 = tolfac_ai
        while(sign(testlb) == sign(testm1)) {
          tolfac_l = tolfac_l + 1
          tolfac_ai1 = tolfac_ai1 + 1
          m1 = round(mu_a_hat - 10 ^ (-tolfac_ai1) * (mu_a_hat - min_a), m1m2lbub_digit)
          lb = min_a + 10 ^ (-tolfac_l) * (mu_a_hat - min_a)
          possibleError4 <- tryCatch(
            testlb <- neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind]),
            error = function(e) e,
            warning=function(w) w
          )
          possibleError4_2 <- tryCatch(
            testm1 <- neg2logR_CB(m1, Tai, EL_CBcrit[as_right_ind]),
            error = function(e) e,
            warning=function(w) w
          )
          if (inherits(possibleError4, "error") | inherits(possibleError4, "warning")) {
            EL_CB[1, ai, as_right_ind] = min_a
            EL_CB_aicheck[1, ai, as_right_ind] = 1
            prev2 = NA
            prev1 = NA
            break
          }
        }
        if (EL_CB[1, ai, as_right_ind] == -1) {
          EL_CB[1, ai, as_right_ind] = uniroot(neg2logR_CB, interval = c(lb, m1), tol = 10 ^ (-tolfac_l), Tai = Tai, EL_CBcrit = EL_CBcrit[as_right_ind])$root
          prev1 = EL_CB[1, ai, as_right_ind] * sc + cen
        }
        if (EL_CB[1,ai, as_right_ind] <= lb) {
          EL_CB[1, ai:n_Ta, as_right_ind] = EL_CB[1, ai, as_right_ind]
          EL_CB_aicheck[1, ai, as_right_ind] = 1
          prev1 = NA
        }
        tolfac_ai2 = tolfac_ai
        while(sign(testub) == sign(testm2)) {
          tolfac_u = tolfac_u + 1
          tolfac_ai2 = tolfac_ai2 + 1
          m2 = round(mu_a_hat + 10 ^ (-tolfac_ai2) * (max_a_orig - mu_a_hat), m1m2lbub_digit)
          ub = max_a_orig - 10 ^ (-tolfac_u) * (max_a_orig - mu_a_hat)
          possibleError5 <- tryCatch(
            testub <- neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind]),
            error = function(e) e,
            warning=function(w) w
          )
          possibleError5_2 <- tryCatch(
            testm2 <- neg2logR_CB(m2, Tai, EL_CBcrit[as_right_ind]),
            error = function(e) e,
            warning=function(w) w
          )
          if (inherits(possibleError5, "error") | inherits(possibleError5, "warning")) {
            EL_CB[2, ai, as_right_ind] = max_a
            EL_CB_aicheck[2, ai, as_right_ind] = 1
            prev2 = NA
            prev1 = NA
            break
          }
        }
        if (EL_CB[2, ai, as_right_ind] == -1) {
          EL_CB[2, ai, as_right_ind] = uniroot(neg2logR_CB, interval = c(m2, ub), tol = 10 ^ (-tolfac_u), Tai = Tai, EL_CBcrit = EL_CBcrit[as_right_ind])$root
          prev2 = EL_CB[2, ai, as_right_ind] * sc + cen
        }
        if (EL_CB[2, ai, as_right_ind] <= lb) {
          EL_CB[2, ai:n_Ta, as_right_ind] = EL_CB[2, ai, as_right_ind]
          EL_CB_aicheck[2, ai, as_right_ind] = 1
          prev2 = NA
        }
        tolfac_ls[ai, as_right_ind] = tolfac_l
        tolfac_us[ai, as_right_ind] = tolfac_u
        if (EL_CB[2, ai, as_right_ind] != -1 & EL_CB[2, ai, as_right_ind] - EL_CB[1, ai, as_right_ind] < (10 ^ - m1m2lbub_digit)) {
          max_a = max(Tai)
          min_a = min(Tai)
          tolfac_ai = 4
          tolfac_l = 10
          tolfac_u = 10
          m1 = mu_a_hat - 10 ^ (-tolfac_ai) * (mu_a_hat - min_a)
          m2 = mu_a_hat + 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
          ub = max_a - 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
          lb = min_a + 10 ^ (-tolfac_ai) * (mu_a_hat - min_a)
          while (sum((Tai - m1) > 0) %in% c(0, n_subject) | sum((Tai - m2) > 0) %in% c(0, n_subject) | sum((Tai - lb) > 0) %in% c(0, n_subject) | sum((Tai - ub) > 0) %in% c(0, n_subject)) {
            tolfac_ai = tolfac_ai + 1
            m1 = mu_a_hat - 10 ^ (-tolfac_ai) * (mu_a_hat - min_a)
            m2 = mu_a_hat + 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
            ub = max_a - 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
            lb = min_a + 10 ^ (-tolfac_ai) * (mu_a_hat - min_a)
          }
          while (sum((Tai - m1) > 0) %in% c(0, n_subject) | sum((Tai - m2) > 0) %in% c(0, n_subject) | sum((Tai - lb) > 0) %in% c(0, n_subject) | sum((Tai - ub) > 0) %in% c(0, n_subject)) {
            tolfac_ai = tolfac_ai + 1
            m1 = round(m1_start - 10 ^ (-tolfac_ai) * (m1_start - min_a), m1m2lbub_digit)
            m2 = round(mu_a_hat + 10 ^ (-tolfac_ai) * (max_a - mu_a_hat), m1m2lbub_digit)
            ub = max_a - 10 ^ (-tolfac_ai) * (max_a - mu_a_hat)
            lb = min_a + 10 ^ (-tolfac_ai) * (m1_start - min_a)
          }
          possibleError11 <- tryCatch(
            testm1 <- neg2logR_CB(m1, Tai, EL_CBcrit[as_right_ind], error_out = 1),
            error = function(e) e
          )
          possibleError12 <- tryCatch(
            testm2 <- neg2logR_CB(m2, Tai, EL_CBcrit[as_right_ind], error_out = 1),
            error = function(e) e
          )
          possibleError2 <- tryCatch(
            testub <- neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind], error_out = 1),
            error = function(e) e
          )
          possibleError3 <- tryCatch(
            testlb <- neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind], error_out = 1),
            error = function(e) e
          )
          if((testm1[2] + testm2[2] + testub[2] + testlb[2] == 0) & !inherits(possibleError11, "error") & !inherits(possibleError12, "error") & !inherits(possibleError2, "error") & !inherits(possibleError3, "error")) {
            testm1 = neg2logR_CB(m1, Tai, EL_CBcrit[as_right_ind])
            testm2 = neg2logR_CB(m2, Tai, EL_CBcrit[as_right_ind])
            testub = neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind])
            testlb = neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind])
            while(sign(testlb) == sign(testm1)) {
              tolfac_l = tolfac_l + 1
              lb = min_a + 10 ^ (-tolfac_l) * (mu_a_hat - min_a)
              possibleError4 <- tryCatch(
                testlb <- neg2logR_CB(lb, Tai, EL_CBcrit[as_right_ind]),
                error = function(e) e,
                warning=function(w) w
              )
              if (inherits(possibleError4, "error") | inherits(possibleError4, "warning")) {
                EL_CB[1, ai, as_right_ind] = min_a
                EL_CB_aicheck[1, ai, as_right_ind] = 1
                prev2 = NA
                prev1 = NA
                break
              }
            }
            EL_CB[1, ai, as_right_ind] = uniroot(neg2logR_CB, interval = c(lb, m1), tol = 10 ^ (-tolfac_l), Tai = Tai, EL_CBcrit = EL_CBcrit[as_right_ind])$root
            if (EL_CB[1,ai, as_right_ind] <= lb) {
              EL_CB[1, ai:n_Ta, as_right_ind] = EL_CB[1, ai, as_right_ind]
              EL_CB_aicheck[1, ai, as_right_ind] = 1
              prev1 = NA
            }
            while(sign(testub) == sign(testm2)) {
              tolfac_u = tolfac_u + 1
              ub = max_a - 10 ^ (-tolfac_u) * (max_a - mu_a_hat)
              possibleError5 <- tryCatch(
                testub <- neg2logR_CB(ub, Tai, EL_CBcrit[as_right_ind]),
                error = function(e) e,
                warning=function(w) w
              )
              if (inherits(possibleError5, "error") | inherits(possibleError5, "warning")){
                EL_CB[2, ai, as_right_ind] = max_a
                EL_CB_aicheck[2, ai, as_right_ind] = 1
                break
              }
            }
            EL_CB[2, ai, as_right_ind] = uniroot(neg2logR_CB, interval = c(m2, ub), tol = 10 ^ (-tolfac_u), Tai = Tai, EL_CBcrit = EL_CBcrit[as_right_ind])$root
            if (EL_CB[2, ai, as_right_ind] <= lb) {
              EL_CB[2, ai:n_Ta, as_right_ind] = EL_CB[2, ai, as_right_ind]
              EL_CB_aicheck[2, ai, as_right_ind] = 1
            }
            tolfac_ls[ai, as_right_ind] = tolfac_l
            tolfac_us[ai, as_right_ind] = tolfac_u
          } else {
            neg2logR_CB_decisionstop_error = 1
            break
          }
        }
      } else {
        neg2logR_CB_decisionstop_error = 1
        break
      }
      EL_CB[, ai, as_right_ind] = EL_CB[, ai, as_right_ind] * sc + cen
    }
  }
  EL_CB_eval = array(-1, c(2, n_agrid_eval, right_endpt_len))
  EL_CBrej = 1:right_endpt_len * NA
  EL_CB_Nairrej = 1:right_endpt_len * NA
  EP_CBrej = 1:right_endpt_len * NA
  EP_CB_Nairrej = 1:right_endpt_len * NA
  HW_CBrej = 1:right_endpt_len * NA
  EL_CB_Nair = list()
  EP_CB_Nair = list()
  for (as_right_ind in 1:right_endpt_len) {
    EL_CB_eval[1, , as_right_ind] = (EL_CB[1, , as_right_ind])[b_a_vec]
    EL_CB_eval[2, , as_right_ind] = (EL_CB[2, , as_right_ind])[b_a_vec]
    EP_CB_Nair[[as_right_ind]] = EP_CB_eval[, , as_right_ind]
    EP_CB_Nair[[as_right_ind]][2, -(1:as_right_ind_keep0[as_right_ind])] = EP_CB_eval[2, as_right_ind_keep0[as_right_ind], as_right_ind]
    EP_as_right = EP_CB_Nair[[as_right_ind]][1, as_right_ind_keep0[as_right_ind]]
    EP_CB_Nair[[as_right_ind]][1, -(1:as_right_ind_keep0[as_right_ind])] = EP_as_right * (EP_as_right < 0)
    EL_CB_Nair[[as_right_ind]] = EL_CB_eval[, , as_right_ind]
    EL_CB_Nair[[as_right_ind]][2, -(1:as_right_ind_keep0[as_right_ind])] = EL_CB_eval[2, as_right_ind_keep0[as_right_ind], as_right_ind]
    EL_CB_Nair[[as_right_ind]][1, -(1:as_right_ind_keep0[as_right_ind])] = 0
    if (!is.na(mu_a[1])) {
      EP_CB_Nairrej[as_right_ind] = (sum(((EP_CB_Nair[[as_right_ind]][2, ] >= mu_a) * (EP_CB_Nair[[as_right_ind]][1, ] <= mu_a)) == 0) >= 1)
      EP_CBrej[as_right_ind] = (sum(((EP_CB_eval[2, , as_right_ind] >= mu_a) * (EP_CB_eval[1, , as_right_ind] <= mu_a)) == 0) >= 1)
      HW_CBrej[as_right_ind] = (sum(((HW_CB_eval[2, , as_right_ind] >= mu_a) * (HW_CB_eval[1, , as_right_ind] <= mu_a)) == 0) >= 1)
      EL_CB_Nairrej[as_right_ind] = (sum(((EL_CB_Nair[[as_right_ind]][2, ] >= mu_a) * (EL_CB_Nair[[as_right_ind]][1, ] <= mu_a)) == 0) >= 1)
      EL_CBrej[as_right_ind] = (sum(((EL_CB_eval[2, , as_right_ind] >= mu_a) * (EL_CB_eval[1, , as_right_ind] <= mu_a)) == 0) >= 1)
    }
  }
  return(list(
    EL_CB_Nair = EL_CB_Nair,
    EL_CB_Nairrej = EL_CB_Nairrej,
    EL_CBcrit = EL_CBcrit,
    EP_CB_Nair = EP_CB_Nair,
    EP_CB_Nairrej = EP_CB_Nairrej,
    EP_CBcrit = EP_CBcrit,
    HW_CB_eval = HW_CB_eval,
    HW_CBrej = HW_CBrej,
    HW_CBcrit = HW_CBcrit,
    mu_hat_vec = mu_hat[1, ],
    S2_hat_vec = S2_hat[1, ],
    EL_CB_aicheck = EL_CB_aicheck,
    neg2logR_CB_decisionstop_error = neg2logR_CB_decisionstop_error
  ))
}

