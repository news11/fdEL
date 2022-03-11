#' Modified \code{\link[fdANOVA:fanova.tests]{fanova.tests}} for Fmax and GPF tests
#'
#' Modifying the Fmaxb (bootstrapped version of Fmax test) and GPF tests computed by \code{\link[fdANOVA]{fanova.tests}} (version 0.1.2) to implement the Uno's selection approach in Supplement Section 5.2
#'
#' @param x a matrix of data, where each column is a discretized version of a function and each row corresponds to each design time point.
#' @param group.label a vector, each element being the group label corresponding to each column of \code{x}.
#' @param n_boot number of bootstrap replications used in calibration.
#' @param as_right_ind_keep_as0 the minimal (among the \code{length(unique(group.label))} groups) index in the vector of design time points corresponding to each element of the set \eqn{{\bf Z}} defined in Supplement Section 5.2.
#' @param parameters a list of parameters for (parallel) computing: dir_path2 is the path to store the (parallel) computing files, no_cores is the number of cores to be used (1 if sequential computing), useseed decides whether one wants to set seed within each core, and seed controls the random number generation used in each core.
#'
#' @return a vector (p-value of the GPF test, p-value of the Fmaxb test).
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @source \href{https://CRAN.R-project.org/package=fdANOVA}{R package fdANOVA} (version 0.1.2)
#'
#' @examples
#' n_subject = sapply(1:4, FUN = function (j) {
#' length(Xt[[j]])
#' })
#' Ta_mat = matrix(unlist(lapply(1:4, FUN = function (j) {
#' c(t(OTdata$Ta[[j]]))
#' })), byrow = TRUE, ncol = length(OTdata$as))
#' anova_onefactor_group = as.factor(rep(1:4, times = n_subject))
#' \dontrun{fanovatests_sel(x = t(Ta_mat), group.label = anova_onefactor_group, as_right_ind_keep_as0 = c(436, 500), parameters = list(dir_path2 = "E:\\R_program\\out", no_cores = 25, useseed = TRUE, seed = 1007))}
fanovatests_sel = function(x = NULL, group.label, n_boot = 1000, as_right_ind_keep_as0 = dim(x)[2], parameters = list(dir_path2 = getwd(), no_cores = 25, useseed = TRUE, seed = 1007)){
  nrep_sub =  n_boot / parameters$no_cores
  split_set =  1 : parameters$no_cores
  right_endpt_len = length(as_right_ind_keep_as0)
  group_sub = as.numeric(unique(group.label))
  division00 <- function(x, y) {
    out <- x/y
    out[as.logical((x == 0)*(x == y))] <- 0
    return (out)
  }
  if(any(is.na(group.label))){ stop("argument group.label can not contain NA values") }
  group.label0 = unique(group.label)
  l = length(group.label0)
  n.i = numeric(l)
  for(i in 1:l) n.i[i] = sum(group.label == group.label0[i])
  x = as.matrix(x); n = ncol(x); p = nrow(x)
  if(n != length(group.label)){
    stop("number of observations (number of columns in x) and number of elements
         in vector of group labels (group.label) must be the same")
  } # END if
  mu0 = rowMeans(x) # length = p
  vmu = matrix(0, nrow = l, ncol = p)
  z = matrix(0, nrow = n, ncol = p)
  SSR = 0; SSE = 0
  for(i in 1:l){
    xi = x[, group.label == group.label0[i]]
    mui = rowMeans(xi); vmu[i,] = mui
    zi = t(xi) - as.matrix(rep(1, n.i[i])) %*% mui
    if(i==1){ z[1:n.i[i],] = zi }else{ z[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = zi }
    SSR = SSR + n.i[i]*(mui - mu0)^2
    SSE = SSE + colSums(zi^2)
  }
  SSR_Nair = list()
  SSE_Nair = list()
  for (as_right_ind in 1:right_endpt_len) {
    SSR_Nair[[as_right_ind]] = SSR
    SSE_Nair[[as_right_ind]] = SSE
    SSR_Nair[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
    SSE_Nair[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
  }
  xs = matrix(0, nrow = n, ncol = p)
  for(i in 1:l){ xs[group.label == group.label0[i],] = t(x[, group.label == group.label0[i]]) - as.matrix(rep(1, n.i[i])) %*% vmu[i,] }
  if (parameters$no_cores > 1) {
    cl <- parallel::makeCluster(parameters$no_cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    foreach::foreach(split = split_set,
      .combine = c
    ) %dopar% {
      if(parameters$useseed == TRUE) {
        set.seed(parameters$seed + sum(group_sub) + split, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      }
      statFmaxboot = matrix(0, nrow = nrep_sub, ncol = right_endpt_len)
      statGPFboot = matrix(0, nrow = nrep_sub, ncol = right_endpt_len)
      for(ii in 1:nrep_sub){
        vmuboot = matrix(0, nrow = l, ncol = p)
        zboot = matrix(0, nrow = n, ncol = p)
        xboot = xs[sample(1:n, replace = TRUE),]
        mu0boot = colMeans(xboot)
        SSRboot = 0
        for(i in 1:l){
          xiboot = xboot[group.label == group.label0[i],]
          muiboot = colMeans(xiboot)
          vmuboot[i,] = muiboot
          ziboot = xiboot - as.matrix(rep(1, n.i[i])) %*% muiboot
          if(i==1){
            zboot[1:n.i[i],] = ziboot
          }else{
            zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot
          }
          SSRboot = SSRboot + n.i[i]*(muiboot - mu0boot)^2
        }
        SSEboot = diag(t(zboot) %*% zboot)
        SSRboot_list = list()
        SSEboot_list = list()
        for (as_right_ind in 1:right_endpt_len) {
          SSRboot_list[[as_right_ind]] = SSRboot
          SSEboot_list[[as_right_ind]] = SSEboot
          SSRboot_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
          SSEboot_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
          statFmaxboot[ii, as_right_ind] = max(division00(SSRboot_list[[as_right_ind]], SSEboot_list[[as_right_ind]])*(n-l)/(l-1))
          statGPFboot[ii, as_right_ind] = mean(division00(SSRboot_list[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]], SSEboot_list[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]])*(n-l)/(l-1))
        }
      }
      setwd(parameters$dir_path2)
      out = list(statFmaxboot = statFmaxboot, statGPFboot = statGPFboot)
      save(out,  file = paste("group_sub_", paste(group_sub, collapse = "_"), "_n_boot_", n_boot, "_split_", split, ".Rdata", sep = ""))
    }
    parallel::stopCluster(cl)
  } else {
    if(parameters$useseed == TRUE) {
      set.seed(parameters$seed + sum(group_sub), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    }
    split = 1
    statFmaxboot = matrix(0, nrow = nrep_sub, ncol = right_endpt_len)
    statGPFboot = matrix(0, nrow = nrep_sub, ncol = right_endpt_len)
    for(ii in 1:nrep_sub){
      vmuboot = matrix(0, nrow = l, ncol = p)
      zboot = matrix(0, nrow = n, ncol = p)
      xboot = xs[sample(1:n, replace = TRUE),]
      mu0boot = colMeans(xboot)
      SSRboot = 0
      for(i in 1:l){
        xiboot = xboot[group.label == group.label0[i],]
        muiboot = colMeans(xiboot)
        vmuboot[i,] = muiboot
        ziboot = xiboot - as.matrix(rep(1, n.i[i])) %*% muiboot
        if(i==1){
          zboot[1:n.i[i],] = ziboot
        }else{
          zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot
        }
        SSRboot = SSRboot + n.i[i]*(muiboot - mu0boot)^2
      }
      SSEboot = diag(t(zboot) %*% zboot)
      SSRboot_list = list()
      SSEboot_list = list()
      for (as_right_ind in 1:right_endpt_len) {
        SSRboot_list[[as_right_ind]] = SSRboot
        SSEboot_list[[as_right_ind]] = SSEboot
        SSRboot_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
        SSEboot_list[[as_right_ind]][-(1:as_right_ind_keep_as0[as_right_ind])] = 0
        statFmaxboot[ii, as_right_ind] = max(division00(SSRboot_list[[as_right_ind]], SSEboot_list[[as_right_ind]])*(n-l)/(l-1))
        statGPFboot[ii, as_right_ind] = mean(division00(SSRboot_list[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]], SSEboot_list[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]])*(n-l)/(l-1))
      }
    }
    setwd(parameters$dir_path2)
    out = list(statFmaxboot = statFmaxboot, statGPFboot = statGPFboot)
    save(out,  file = paste("group_sub_", paste(group_sub, collapse = "_"), "_n_boot_", n_boot, "_split_", split, ".Rdata", sep = ""))
  }
  split_ind = 1
  statFmaxboot = matrix(0, nrow = n_boot, ncol = right_endpt_len)
  statGPFboot = matrix(0, nrow = n_boot, ncol = right_endpt_len)
  setwd(parameters$dir_path2)
  for (split in split_set) {
    load(paste("group_sub_", paste(group_sub, collapse = "_"), "_n_boot_", n_boot, "_split_", split, ".Rdata", sep = ""))
    statFmaxboot[(nrep_sub * (split_ind - 1) + (1:nrep_sub)), ] = out$statFmaxboot
    statGPFboot[(nrep_sub * (split_ind - 1) + (1:nrep_sub)), ] = out$statGPFboot
    split_ind = split_ind + 1
  }
  p_c_Fmaxb = 1:right_endpt_len * 0
  P_c_Fmaxb = matrix(0, nrow = right_endpt_len, ncol = n_boot)
  p_c_GPF = 1:right_endpt_len * 0
  P_c_GPF = matrix(0, nrow = right_endpt_len, ncol = n_boot)
  for (as_right_ind in 1:right_endpt_len) {
    statFmax = max(division00(SSR_Nair[[as_right_ind]], SSE_Nair[[as_right_ind]])*(n-l)/(l-1))
    statGPF = mean(division00(SSR_Nair[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]], SSE_Nair[[as_right_ind]][1:as_right_ind_keep_as0[as_right_ind]])*(n-l)/(l-1))
    p_c_Fmaxb[as_right_ind]=mean(statFmaxboot[, as_right_ind]>statFmax) # vector of length(c_seq)
    P_c_Fmaxb[as_right_ind, ]=(n_boot-rank(statFmaxboot[, as_right_ind],ties.method ="max"))/n_boot
    p_c_GPF[as_right_ind]=mean(statGPFboot[, as_right_ind]>statGPF) # vector of length(c_seq)
    P_c_GPF[as_right_ind, ]=(n_boot-rank(statGPFboot[, as_right_ind],ties.method ="max"))/n_boot
  }
  p_b_Fmaxb=min(p_c_Fmaxb)
  P_b_Fmaxb=apply(P_c_Fmaxb,2,min)
  pvalueFmaxb = mean(p_b_Fmaxb > P_b_Fmaxb)
  p_b_GPF=min(p_c_GPF)
  P_b_GPF=apply(P_c_GPF,2,min)
  pvalueGPF = mean(p_b_GPF > P_b_GPF)
  return(c(pvalueGPF, pvalueFmaxb))
}
