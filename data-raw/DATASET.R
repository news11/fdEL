## read in raw data:
rawact <- readRDS("data-raw/rawact_load.Rdata")
demog <- readRDS("data-raw/demog_load.Rdata")

## prepare OTdata and rawlens (occupation time related) dataset:
attach(demog)
band_select_groups = list()
band_select_groups[[1]] = seqn[(!is.na(ridageyr) & ridageyr >= 75 & ridageyr <= 85 & !is.na(dmqmilit) & dmqmilit == 1)]
band_select_groups[[2]] = seqn[(!is.na(ridageyr) & ridageyr >= 75 & ridageyr <= 85 & !is.na(dmqmilit) & dmqmilit == 2)]
band_select_groups[[3]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 & ridageyr < 75 & !is.na(dmqmilit) & dmqmilit == 1)]
band_select_groups[[4]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 & ridageyr < 75 & !is.na(dmqmilit) & dmqmilit == 2)]

n_group = length(band_select_groups)
n_subject = 1:n_group * 0  # ok
rawact_band_select = list()

Xt = sapply(1:n_group, FUN = function (j) {
  rawact_band_select[[j]] <<- rawact[rawact[, 1] %in% band_select_groups[[j]], ]
  attach(rawact_band_select[[j]])
  subject_ids = unique(seqn)
  n_subject[j] <<- length(unique(seqn))
  t_ncol = length(unique(paxn))
  lapply(1:n_subject[j], FUN = function (i) {
    subject_id = subject_ids[i]
    paxinten[seqn == subject_id & paxstat==1 & paxcal==1]
  })  # paxinten is the activity count; see, e.g., the paper ``Organizing and analyzing the activity data in NHANES'' by Leroux et al.
})
usethis::use_data(Xt, overwrite = TRUE)

rawlens = list()
as = seq(0, 499, by = parameters$by_agrid)
mu_a_hat = function(Xt, grid_widths, as) {
  # Computes the activity profile of a stochastic process of accelerometer readings
  # Args:
  #   Xt: a vector of accelerometer readings on a grid of time points
  #   grid_widths: the number of time units each grid point contains
  #   Xt, grid_widths should be vectors of the same length
  #   as: the grid of a's
  # Returns:
  #   activity profile of Xt
  n_a = length(as)
  act_profile = sapply(1:n_a, FUN = function (a) {
    sum(grid_widths[Xt > as[a]])
  })
  return(act_profile)
}
Ta = lapply(1:n_group, FUN = function (j) {
  observed_Xt_lens_tmp = 1:n_subject[j] * 0
  Ta_j = t(sapply(1:n_subject[j], FUN = function (i) {
    observed_Xt = (Xt[[j]])[[i]][!is.na((Xt[[j]])[[i]])]
    observed_Xt_len = length(observed_Xt)
    observed_Xt_lens_tmp[i] <<- observed_Xt_len
    grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)
    mu_a_hat(observed_Xt, grid_widths = grid_width_i, as)
  }))
  rawlens[[j]] <<- observed_Xt_lens_tmp
  return(Ta_j)
})

OTdata = list(Ta = Ta, as = as)
usethis::use_data(OTdata, overwrite = TRUE)
usethis::use_data(rawlens, overwrite = TRUE)

############ to-do: simulated data
# parameters:
# [0, tau]: the time domain of Xt process
parameters$tau = 1
parameters$n_floor = 1
parameters$by_agrid = 1 / parameters$n_floor * 1
floor_n = function(x, n = parameters$n_floor) {
  # A more flexible floor function.
  # The round function is to eliminate a small numerical difference from the actual value of x * n.
  # Args:
  #   x: a number such that m / n <= x <= (m + 1) / n for a non-negative interger m
  # Returns:
  #   m / n
  floor(round(x * n, 2)) / n
}

####### data too large that slows the pkg down, so give up
## prepare (raw accelerometer) dataset:
demog_use <- as.data.frame(cbind(demog$seqn, demog$ridageyr, demog$dmqmilit))

usethis::use_data(demog_use, overwrite = TRUE)
rawact_use <- as.data.frame(cbind(rawact$seqn, rawact$paxstat, rawact$paxcal, rawact$paxn, rawact$paxinten  ))

usethis::use_data(rawact_use, overwrite = TRUE)

