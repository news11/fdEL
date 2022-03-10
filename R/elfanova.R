#' EL and related functional ANOVA tests
#'
#' Computes EL-based functional ANOVA tests and the corresponding asymptotically equivalent Wald-type tests
#'
#' @param Ta a list with the j-th element being the matrix of observed functional data (e.g., occuptation time curve) for the j-th group, and the i-th row of that matrix representing the observed curve of the i-th subject in the j-th group
#' @param as a vector consisting of the design time points (e.g., activity levels) of the observed functional data; that is, \eqn{{\bf G}_n} defined in Section 2.2 (from the smallest to the largest).
#' @param n_boot number of bootstrap samples.
#' @param as_right_ind_keep_as0 the minimal (among the \code{length(Ta)} groups) index in the vector \code{as} (defined above) corresponding to each element of the set \eqn{M} (defined in Suppplement Section 5.2).
#'
#' @return a list containing the following:
#'   \item{suptest}{test statistic of the proposed maximally selected EL test}
#'   \item{inttest_da}{test statistic of the integrated EL test}
#'   \item{suptest_EP}{test statistic of the Wald-type test that is asymptotically equivalent to the proposed maximally selected EL test}
#'   \item{inttest_EP_da}{test statistic of the Wald-type test that is asymptotically equivalent to the integrated EL test}
#'   \item{out_sup_pval}{p-value of the proposed maximally selected EL test}
#'   \item{out_dF_pval}{p-value of the integrated EL test}
#'   \item{out_sup_EP_pval}{p-value of the Wald-type test that is asymptotically equivalent to the proposed maximally selected EL test}
#'   \item{out_dF_EP_pval}{p-value of the Wald-type test that is asymptotically equivalent to the integrated EL test}
#'   \item{as}{the input \code{as}; see Arguments}
#'   \item{mu_hat_vec}{a \code{length(Ta)} x \code{length(as)} matrix; j-th row is the sample mean of \code{Ta[[j]]}}
#'   \item{supboot}{bootstrap values of the proposed maximally selected EL test}
#'   \item{int_da_boot}{bootstrap values of the integrated EL test}
#'
#' @export
#'
#' @examples
#' set.seed(1008, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
#' test_out = elfanova(Ta = OTdata$Ta, as = OTdata$as)
#' test_out$out_sup_pval
elfanova = function(Ta, as, n_boot = 1000, as_right_ind_keep_as0 = length(as)) {
  n_Ta = length(as) # the length of as
  n_subject = sapply(Ta, FUN = function (x) {
    dim(x)[1]
  }) # a vector, n_subject[j] being the sample size for the j-th group
  gamma_js = n_subject / sum(n_subject)
  U_hat_star = array(0, c(n_boot, n_Ta, length(n_subject)))
  theta_hat_js = array(0, c(n_boot, n_Ta, length(n_subject)))
  wjs = array(0, c(n_boot, n_Ta, length(n_subject)))
  mu_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  S2_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  boot_indx = list()
  U_hat_star_run = sapply(1:length(n_subject), FUN = function (j) {
    boot_indx[[j]] <<- sample(1:n_subject[j], n_subject[j] * n_boot, replace = T) # suppose n_subject[j] subjects selected, then repeated n_boot times
    boot_indx_mat = matrix(boot_indx[[j]], byrow = T, nrow = n_boot, ncol = n_subject[j])
    mu_hat = matrix(rep(apply(Ta[[j]], 2, mean), each = n_boot), nrow = n_boot, ncol = n_Ta)
    mu_hat_vec[j, ] <<- mu_hat[1, ]
    Ta_boot = array((Ta[[j]])[boot_indx[[j]], ], c(n_subject[j], n_boot, n_Ta))
    mu_hat_star = apply(Ta_boot, c(2, 3) , mean)
    S2_hat = matrix(apply(Ta[[j]], 2, var) * (n_subject[j] - 1) / n_subject[j], byrow = T, nrow = n_boot, ncol = n_Ta)
    S2_hat_vec[j, ] <<- S2_hat[1, ]
    U_hat_star[, , j] <<- sqrt(n_subject[j]) * U_division00_mat(mu_hat_star -  mu_hat, sqrt(S2_hat)) # n_boot x n_Ta
    theta_hat_js[, , j] <<- S2_hat / gamma_js[j]
    wjs[, , j] <<- 1 / theta_hat_js[, , j] # need to be normalized to 1 later; but just once for all bootstrap samples
    return(0) # like Ujps in ELSOc_k
  })
  wjs = product_arr3(wjs, 1 / array(rep(apply(wjs[1, , ], 1, sum), each = n_boot), c(n_boot, n_Ta, length(n_subject))))
  wjs_vec = t(wjs[1, , ])
  avg_U_hat_star = apply(sqrt(wjs) * U_hat_star, c(1, 2), sum)
  avg_U_hat_star_karray = array(rep(avg_U_hat_star, times = length(n_subject)), c(n_boot, n_Ta, length(n_subject)))
  U2_boot = apply(wjs * (product_arr3(U_hat_star, 1 / sqrt(wjs)) - avg_U_hat_star_karray) ^ 2, c(1, 2), sum) # bootstrap SSB(a) in Theorem 1, n_boot x n_Ta

### EL test:
  teststat_pre = 1:n_Ta * 0
  error_vec = 1:n_Ta * 0
  for (ai in 1:max(as_right_ind_keep_as0)) {

    neg2logRa = neg2logR_test(ai, Ta, n_subject, EL_testcrit = 0, error_out = 1)
    error_vec[ai] = neg2logRa$still_error
    if (error_vec[ai] == 1) next
    teststat_pre[ai]=neg2logRa$neg2logR_crit
  }

# Wald-type:
  mu_bar = apply(matrix(gamma_js, nrow = length(n_subject), ncol = n_Ta) * mu_hat_vec, 2, sum)
  mu_bar_vec = matrix(mu_bar, byrow = T, nrow = length(n_subject), ncol = n_Ta)
  n_subject_vec = matrix(n_subject, nrow = length(n_subject), ncol = n_Ta)
  U_hat = sqrt(n_subject_vec) * U_division00_mat(mu_hat_vec -  mu_bar_vec, sqrt(S2_hat_vec))
  avg_U_hat = apply(sqrt(wjs_vec) * U_hat, 2, sum)
  avg_U_hat_vec = matrix(avg_U_hat, byrow = T, nrow = length(n_subject), ncol = n_Ta)
  U2_vec = apply(wjs_vec * (U_division00_mat(U_hat, sqrt(wjs_vec)) - avg_U_hat_vec) ^ 2, 2, sum) # estimate of SSB(a) in Theorem 1, n_boot x n_Ta

  right_endpt_len = length(as_right_ind_keep_as0)
  EL_Nair_sup = list()
  EL_Nair_int = list()
  EP_Nair_sup = list()
  EP_Nair_int = list()
  sup_boot = list()
  int_da_boot = list()
  p_c_EP_Nair_sup = 1:right_endpt_len * 0
  p_c_EP_Nair_int = 1:right_endpt_len * 0
  p_c_EL_Nair_sup = 1:right_endpt_len * 0
  p_c_EL_Nair_int = 1:right_endpt_len * 0
  U2_vec_list = list()
  teststat_pre_list = list()
  U2_boot_list = list()
  P_c_sup = matrix(0, nrow = right_endpt_len, ncol = n_boot)
  P_c_int = matrix(0, nrow = right_endpt_len, ncol = n_boot)
  for (as_right_ind in 1:right_endpt_len) {
    U2_vec_list[[as_right_ind]] = U2_vec
    U2_vec_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
    teststat_pre_list[[as_right_ind]] = teststat_pre
    teststat_pre_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
    EP_Nair_sup[[as_right_ind]] = max(U2_vec_list[[as_right_ind]])
    EP_Nair_int_pre_da = U2_vec_list[[as_right_ind]][-n_Ta] * diff(as)
    EP_Nair_int[[as_right_ind]] = sum(EP_Nair_int_pre_da)
    EL_Nair_sup[[as_right_ind]] = max(teststat_pre_list[[as_right_ind]])
    EL_Nair_int_pre_da = teststat_pre_list[[as_right_ind]][-n_Ta] * diff(as)
    EL_Nair_int[[as_right_ind]] = sum(EL_Nair_int_pre_da)
    U2_boot_list[[as_right_ind]] = U2_boot
    U2_boot_list[[as_right_ind]][, -(1:as_right_ind_keep_as0[as_right_ind])] = 0
    sup_boot[[as_right_ind]] = apply(as.matrix(U2_boot_list[[as_right_ind]]), 1, max)
    a_big = matrix(rep(as, times = n_boot), byrow = TRUE, nrow = n_boot) # n_boot x n_Ta
    U2_boot_times_da = as.matrix(U2_boot_list[[as_right_ind]][, -n_Ta]) * as.matrix(t(apply(a_big, 1, diff)))
    int_da_boot[[as_right_ind]] = apply(as.matrix(U2_boot_times_da), 1, sum)
    p_c_EP_Nair_sup[as_right_ind]=mean(sup_boot[[as_right_ind]]>EP_Nair_sup[[as_right_ind]]) # vector of length(c_seq)
    p_c_EP_Nair_int[as_right_ind]=mean(int_da_boot[[as_right_ind]]>EP_Nair_int[[as_right_ind]]) # vector of length(c_seq)
    p_c_EL_Nair_sup[as_right_ind]=mean(sup_boot[[as_right_ind]]>EL_Nair_sup[[as_right_ind]]) # vector of length(c_seq)
    p_c_EL_Nair_int[as_right_ind]=mean(int_da_boot[[as_right_ind]]>EL_Nair_int[[as_right_ind]]) # vector of length(c_seq)
    P_c_sup[as_right_ind, ]=(n_boot-rank(sup_boot[[as_right_ind]],ties.method ="max"))/n_boot
    P_c_int[as_right_ind, ]=(n_boot-rank(int_da_boot[[as_right_ind]],ties.method ="max"))/n_boot
  }
  p_b_EP_Nair_sup=min(p_c_EP_Nair_sup)
  p_b_EP_Nair_int=min(p_c_EP_Nair_int)
  p_b_EL_Nair_sup=min(p_c_EL_Nair_sup)
  p_b_EL_Nair_int=min(p_c_EL_Nair_int)
  P_b_sup=apply(P_c_sup,2,min)
  P_b_int=apply(P_c_int,2,min)
  out_sup_pval = mean(p_b_EL_Nair_sup > P_b_sup)
  out_dF_pval = mean(p_b_EL_Nair_int > P_b_int)
  out_sup_EP_pval = mean(p_b_EP_Nair_sup  > P_b_sup)
  out_dF_EP_pval = mean(p_b_EP_Nair_int > P_b_int)

  return(list(
    suptest = EL_Nair_sup,
    inttest_da = EL_Nair_int,
    suptest_EP = EP_Nair_sup,
    inttest_EP_da = EP_Nair_int,
    out_sup_pval = out_sup_pval,
    out_dF_pval = out_dF_pval,
    out_sup_EP_pval = out_sup_EP_pval,
    out_dF_EP_pval = out_dF_EP_pval,
    as = as,
    mu_hat_vec = mu_hat_vec,
    sup_boot = sup_boot,
    int_da_boot = int_da_boot
  ))
}

