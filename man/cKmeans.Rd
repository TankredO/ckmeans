% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckmeans.R
\name{ckmeans}
\alias{ckmeans}
\title{Consensus K-means clustering}
\usage{
ckmeans(
  x,
  k,
  n_rep = 50,
  p_pred = 1,
  p_samp = 1,
  save_kms = TRUE,
  hclust_options = list(method = "average"),
  calc_bic = TRUE,
  ...
)
}
\arguments{
\item{x}{matrix with samples in rows and features in columns}

\item{k}{number of clusters}

\item{n_rep}{number of individual k means runs}

\item{p_pred}{proportion of predictors used in every k means run}

\item{p_samp}{proportion of samples used in every k means run}

\item{save_kms}{logical, (or 'minimize') determining whether the k means object should be saved. This can be very memory demanding, depending on n_rep and x. If 'minimize', the kmeans objects are saved without column names, saving memory.}

\item{hclust_options}{list of option passed to hclust, which is used to generate consensus clusters}

\item{calc_bic}{logical, determining whether the BIC (Bayesian Information Criterion) should be calculated for the k means runs}

\item{...}{arguments passed to kmeans}
}
\value{
ckmeans object
}
\description{
Consensus k means clustering
}
\details{
Runs several independent k means clustering steps, and combines information from the different runs to calculate consensus clusters using hierarchical clustering. The hierarchical clustering is based on the proportion of runs in which each pair of samples is placed in the same cluster, interpreted as distance.
}
\examples{
## generate data
x1 = c(rnorm(10), rnorm(14, 5, 1) + 2, rnorm(30, -5, 1))
x2 = c(rnorm(10), rnorm(14, 5, 1) + 2, rnorm(30, -5, 1))
x3 = c(rnorm(10)-3, rnorm(14, 2, 1) + 2, rnorm(30, 1, 1))
x = matrix(c(x1, x2, x3), ncol = 3, dimnames = list(1:54, c('x1', 'x2', 'x3')))

pairs(x)

## run ckmeans for a single K
ckm = ckmeans(x, 3, n_rep = 100, p_samp = 0.5, p_pred = 0.5)

# plot consensus matrix with color coded clusters
plot(ckm, cex.axis = 0.75)

plotDist(ckm)



# plot(x, col=ckm$cc, pch=c(rep(1, 10), rep(2, 14)))


## run ckmeans for multiple K
ckms = multickmeans(x, 1:7, n_rep = 100, p_samp = 0.8, p_pred = 0.5)
plot(ckms$bics, type='l')
plot(ckms$aics, type='l')
plot(ckms$sils, type='l')
plot(ckms$dbs, type='l')

ckms$

ckms$aics

for (i in 1:length(ckms$ckms)) {
  plotDist(ckms$ckms[[i]], ord=TRUE)
}
}
\author{
Tankred Ott
}
