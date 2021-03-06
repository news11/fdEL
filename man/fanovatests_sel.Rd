% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fanovatests_sel.R
\name{fanovatests_sel}
\alias{fanovatests_sel}
\title{Modified \code{\link[fdANOVA:fanova.tests]{fanova.tests}} for Fmax and GPF tests}
\source{
\href{https://CRAN.R-project.org/package=fdANOVA}{R package fdANOVA} (version 0.1.2)
}
\usage{
fanovatests_sel(
  x = NULL,
  group.label,
  n_boot = 1000,
  as_right_ind_keep_as0 = dim(x)[2],
  parameters = list(dir_path2 = getwd(), no_cores = 25, useseed = TRUE, seed = 1007)
)
}
\arguments{
\item{x}{a matrix of data, where each column is a discretized version of a function and each row corresponds to each design time point.}

\item{group.label}{a vector, each element being the group label corresponding to each column of \code{x}.}

\item{n_boot}{number of bootstrap replications used in calibration.}

\item{as_right_ind_keep_as0}{the minimal (among the \code{length(unique(group.label))} groups) index in the vector of design time points corresponding to each element of the set \eqn{{\bf Z}} defined in Supplement Section 5.2.}

\item{parameters}{a list of parameters for (parallel) computing: dir_path2 is the path to store the (parallel) computing files, no_cores is the number of cores to be used (1 if sequential computing), useseed decides whether one wants to set seed within each core, and seed controls the random number generation used in each core.}
}
\value{
a vector (p-value of the GPF test, p-value of the Fmaxb test).
}
\description{
Modifying the Fmaxb (bootstrapped version of Fmax test) and GPF tests computed by \code{\link[fdANOVA]{fanova.tests}} (version 0.1.2) to implement the Uno's selection approach in Supplement Section 5.2
}
\examples{
n_subject = sapply(1:4, FUN = function (j) {
length(Xt[[j]])
})
Ta_mat = matrix(unlist(lapply(1:4, FUN = function (j) {
c(t(OTdata$Ta[[j]]))
})), byrow = TRUE, ncol = length(OTdata$as))
anova_onefactor_group = as.factor(rep(1:4, times = n_subject))
\dontrun{fanovatests_sel(x = t(Ta_mat), group.label = anova_onefactor_group, as_right_ind_keep_as0 = c(436, 500), parameters = list(dir_path2 = "E:\\\\R_program\\\\out", no_cores = 25, useseed = TRUE, seed = 1007))}
}
