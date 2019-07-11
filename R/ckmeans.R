#' Consensus k means clustering
#' @title Consensus K-means clustering
#' @author Tankred Ott
#' @details Runs several independent k means clustering steps, and combines information from the different runs to calculate consensus clusters using hierarchical clustering. The hierarchical clustering is based on the proportion of runs in which each pair of samples is placed in the same cluster, interpreted as distance.
#' @param x matrix with samples in rows and features in columns
#' @param k number of clusters
#' @param n_rep number of individual k means runs
#' @param p_pred proportion of predictors used in every k means run
#' @param p_samp proportion of samples used in every k means run
#' @param save_kms logical, (or 'minimize') determining whether the k means object should be saved. This can be very memory demanding, depending on n_rep and x. If 'minimize', the kmeans objects are saved without column names, saving memory.
#' @param hclust_options list of option passed to hclust, which is used to generate consensus clusters
#' @param calc_bic logical, determining whether the BIC (Bayesian Information Criterion) should be calculated for the k means runs
#' @param ... arguments passed to kmeans
#' @return ckmeans object
#' @example /inst/examples/ckmeans_example.R
#' @import stats
#' @export
ckmeans <- function(x, k, n_rep = 50, p_pred = 1, p_samp = 1, save_kms = TRUE, hclust_options = list(method = 'average'),
                    calc_bic = TRUE, ...) {

  pred_sel <- NULL
  samp_sel <- NULL
  sampled_together <- NULL
  kms <- NULL
  bics <- NULL

  minimize <- if (save_kms %in% c('minimize', 'min', 'm')) TRUE else FALSE
  if (minimize) save_kms <- TRUE

  if (is.null(hclust_options$method)) hclust_options$method <- 'complete'

  par_names <- colnames(x)
  samp_names <- row.names(x)

  # prepare data structures
  if (p_pred != 1) {
    n_pred <- floor(ncol(x) * p_pred)
    pred_sel <- matrix(0, ncol = n_pred, nrow = n_rep) # to store selected predictors at each rep
  }

  if (p_samp != 1) {
    n_samp <- floor(nrow(x) * p_pred)
    sampled_together <- matrix(0, ncol = nrow(x), nrow = nrow(x)) # to count times x1,...,xn are sampled together
    samp_sel <- matrix(0, ncol = n_samp, nrow = n_rep) # to store selected samples at each rep
  }

  if (save_kms != FALSE) {
    kms <- vector(mode = 'list', length = n_rep)
    if (calc_bic == TRUE) bics <- vector(mode = 'numeric', length = n_rep)
  }
  clustered_together <- matrix(0, ncol = nrow(x), nrow = nrow(x)) # to count times x1,...,xn are clustered together

  # run n_rep kmeans
  for(i in 1:n_rep) {
    .pred_sel <- if (p_pred == 1) 1:ncol(x) else sample(ncol(x), n_pred) # select predictors
    .samp_sel <- if (p_samp == 1) 1:nrow(x) else sample(nrow(x), n_samp) # select samples

    if (p_pred != 1){
      pred_sel[i,] <- .pred_sel # store selected predictors
    }

    if (p_samp != 1){
      sampled_together[.samp_sel, .samp_sel] <- sampled_together[.samp_sel, .samp_sel] + 1
      samp_sel[i,] <- .samp_sel # store selected samples
    }

    # k means
    .km <- kmeans(x = x[.samp_sel, .pred_sel], centers = k, ...)
    .cl <- .km$cluster

    # store number of times samples are clustered together
    if (p_samp == 1) clustered_together <- clustered_together + sapply(.cl, function(x) as.numeric(x == .cl))
    else clustered_together[samp_sel[i,], samp_sel[i,]] <-
      clustered_together[samp_sel[i,], samp_sel[i,]] + sapply(.cl, function(x) as.numeric(x == .cl))

    # save kmeans object
    if (save_kms != FALSE) {
      if (calc_bic == TRUE) bics[i] <- bic_kmeans(.km)
      if (minimize) {
        names(.km$cluster) <- NULL
        colnames(.km$centers) <- NULL
      }
      kms[[i]] <- .km
    }
  }

  # calculate proportion of shared clusters for every pair of individuals
  pcc <- if (p_samp == 1) clustered_together / n_rep else clustered_together / sampled_together
  colnames(pcc) <- samp_names

  # calculate consensus clusters
  d <- as.dist(1-pcc)
  hclust_options$d <- d
  hc <- do.call(hclust, args = hclust_options)
  cc <- cutree(hc, k)

  out <- list(pcc = pcc,
              cc = cc,
              kms = kms,
              bics = bics,
              k = k,
              n_rep = n_rep,
              p_pred = p_pred,
              p_samp = p_samp,
              pred_sel = pred_sel,
              samp_sel = samp_sel,
              sampled_together = sampled_together,
              clustered_together = clustered_together,
              par_names = par_names,
              samp_names = samp_names,
              x = x)
  class(out) <- c('ckmeans')
  out$bic <- bic_kmeans2(out$x, out$cc)
  out
}

#' Consensus k means clustering for multiple ks
#' @title Consensus K-means clustering over mutiple K
#' @author Tankred Ott
#' @details Runs several independent k means clustering steps for several K, and combines information from the different runs to calculate consensus clusters using hierarchical clustering. The hierarchical clustering is based on the proportion of runs in which each pair of samples is placed in the same cluster, interpreted as distance.
#' @param x matrix with samples in rows and features in columns
#' @param ks vector of cluster numbers
#' @param n_rep number of individual k means runs
#' @param p_pred proportion of predictors used in every k means run
#' @param p_samp proportion of samples used in every k means run
#' @param save_kms logical, (or 'minimize') determining whether the k means object should be saved. This can be very memory demanding, depending on n_rep and x. If 'minimize', the kmeans objects are saved without column names, saving memory.
#' @param hclust_options list of option passed to hclust, which is used to generate consensus clusters
#' @param calc_bic logical, determining whether the BIC (Bayesian Information Criterion) should be calculated for the k means runs
#' @param ... arguments passed to kmeans
#' @return multickmeans object. ckms contains the list of ckmeans objects, ks is the vector of ks, bics is the vector of BICs.
#' @export
multickmeans <- function(x, ks, n_rep = 50, p_pred = 1, p_samp = 1, save_kms = TRUE, hclust_options = list(),
                         calc_bic = TRUE, ...) {
  # run consensus k means
  ckms <- vector('list', length(ks))
  for(i in 1:length(ks)) {
    ckms[[i]] <- ckmeans(x = x, k = ks[i], n_rep = n_rep, p_pred = p_pred, p_samp = p_samp, save_kms = save_kms,
                         hclust_options = hclust_options, calc_bic = calc_bic, ...)
  }

  # calculate BIC for consensus clustering
  bics <- unlist(lapply(ckms, function(ckm) ckm$bic))
  names(bics) <- ks

  list(
    ckms = ckms,
    ks = ks,
    bics = bics
  )
}


## Helper functions

# https://en.wikipedia.org/wiki/Bayesian_information_criterion: "Gaussian special case"
#' Function to calculate the BIC (Bayesian Information Criterion) for a kmeans object
#' @param fit kmeans object returned from kmeans()
bic_kmeans <- function(fit){
  n <- length(fit$cluster)
  k <- nrow(fit$centers)
  D <- fit$tot.withinss
  return(n * log(D/n) + log(n)*k)
}

# adapted from adegenet
#' Function to calculate the BIC (Bayesian Information Criterion) for a matrix and cluster membership
#' @param x matrix with samples as rows and predictors as columns
#' @param cl vector of cluster memberships for the samples
bic_kmeans2 <- function (x, cl) {
  cs <- unique(cl)
  k <- length(unique(cl))
  N <- nrow(x)
  RSS <- if (k == 1) sum(apply(x, 2, function(v) sum((v - mean(v))^2))) else .compute.wss(x, cl)

  bic <- N * log(RSS/N) + log(N) * k
  bic
}

#' Helper function for bic_kmeans2
#' @param x matrix with samples as rows and predictors as columns
#' @param f vector of cluster memberships for the samples
.compute.wss <- function(x, f) {
  x.group.mean <- apply(x, 2, tapply, f, mean)
  sum((x - x.group.mean[as.character(f),])^2)
}
