% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elfband.R
\name{elfband}
\alias{elfband}
\title{EL and related confidence bands for functional means}
\usage{
elfband(
  Ta,
  as,
  as_eval,
  n_boot = 1000,
  mu_a = NA,
  alpha = 0.05,
  as_right_ind_keep0 = length(as_eval),
  as_right_ind_keep_as0 = length(as)
)
}
\arguments{
\item{Ta}{a matrix with the i-th row representing the i-th subject's observed functional data (e.g., occupation time data), the j-th column corresponding to the observed functional data at the j-th design time point (e.g., activity level) for each subject}

\item{as}{a vector consisting of the design time points of the observed functional data; that is, \eqn{{\bf G}_n} defined in Section 2.2 (from the smallest to the largest).}

\item{as_eval}{a vector consisting of the time points (from the smallest to the largest; e.g., the activity levels of interest) at which we want to evaluate the confidence band. This vector can be denser or coarser than \code{as}, as long as the largest point of \code{as_eval} is no greater than the last point of \code{as} (defined above).}

\item{n_boot}{number of bootstrap samples.}

\item{mu_a}{a vector containing the true functional mean over \code{as_eval} (defined above)}

\item{alpha}{a number, the significance level of interest}

\item{as_right_ind_keep0}{This is the index in the vector \code{as_eval} (defined above) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1. The default value is \code{length(as_eval)}, where the approach in Supplement Section 5.1 is NOT implemented.}

\item{as_right_ind_keep_as0}{This is the index in the vector \code{as} (defined above) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1. The default value is \code{length(as)}, where the approach in Supplement Section 5.1 is NOT implemented.}
}
\value{
a list containing the following:
\item{EL_CB_Nair}{EL confidence band over \code{as_eval} (defined above) after adapting Nair's two-step approach described in Supplement Section 5.1, with first row being the lower bound and second row being the upper bound}
\item{EL_CB_Nairrej}{1 if the truth \code{mu_a} is not included in the EL confidence band; 0 otherwise}
\item{EL_CBcrit}{critical value for computing the EL confidence band}
\item{EP_CB_Nair}{EP confidence band over \code{as_eval} (defined above) after adapting Nair's two-step approach described in Supplement Section 5.1, with first row being the lower bound and second row being the upper bound}
\item{EP_CB_Nairrej}{1 if the truth \code{mu_a} is not included in the EP confidence band; 0 otherwise}
\item{EP_CBcrit}{critical value for computing the EP confidence band}
\item{HW_CB_eval}{NS confidence band over \code{as_eval} (defined above), with first row being the lower bound and second row being the upper bound}
\item{HW_CBrej}{1 if the truth \code{mu_a} is not included in the NS confidence band; 0 otherwise}
\item{HW_CBcrit}{critical value for computing the NS confidence band}
\item{mu_hat_vec}{a \code{length(as)}-vector containing \eqn{\hat{\mu}(a)}, the sample mean of the observed functional data}
\item{S2_hat_vec}{a \code{length(as)}-vector containing \eqn{\hat{S}^2(a)}, the estimated variance of the observed functional data defined in Section 2.3}
\item{EL_CB_aicheck}{an array, whose (bi, ai, zi) entry recording the following information for the bi-th bound (bi = 1 for lower bound, bi = 2 for upper bound) of the EL band at \code{as[ai]} when adapting Nair's two-step approach described in Supplement Section 5.1 using the zi-th element in the vector of \code{as_right_ind_keep_as0} (defined above): the entry is 1 if the bound EL band is no greater than the minimum of the observed functional data at \code{as[ai]} plus some tolerance value, or when the bound of the EL band is too close to the extreme values of the observed functional data at \code{as[ai]}; 0 otherwise}
\item{neg2logR_CB_decisionstop_error}{1 if there is an error during the EL confidence band construction; 0 otherwise}
}
\description{
Computes EL-based confidence bands and related Wald-type confidence bands for the mean of the functional data of interest, along with whether they capture the true functional mean or not
}
\examples{
CB_out <- elfband(Ta = OTdata$Ta[[1]], as = OTdata$as, as_eval = OTdata$as)
CB_out$EL_CB_Nair
}
