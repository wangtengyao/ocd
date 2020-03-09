#' Processing a new data point for the 'OCD' class
#' @description This function implements the \code{\link{getData}} function to
#' perform the online changepoint detection for the 'OCD' class.
#' @param x_new a new data point
#' @param A matrix of tail CUSUMs that will be tracked and updated online
#' @param tail matrix of tail lengths that will be tracked and updated online
#' @param beta a user specified lower bound on the l_2 norm of the vector of
#' change to be detected.
#' @param sparsity a user specified mode parameter. If the vector of change is
#' known to be dense or sparse, then one should set sparsity to 'dense' or 'sparse'
#' accordingly, otherwise, the default choice sparsity='auto' will run the algorithm
#' adaptive to the sparsity level.
#' @return a list of
#' \itemize{
#' \item stat: a vector of the test statistics for the 'OCD' class
#' \item A: the updated A matrix
#' \item tail: the updated tail matrix
#' }
#' @export
ocd_update <- function(x_new, A, tail, beta, sparsity){
  p <- length(x_new)   # dimension of data
  a <- sqrt(2*log(p))  # hard threshold parameter
  L <- floor(log2(p))  # number of different scales
  tmp <- beta/sqrt(2^(0:(L+1)) * log2(2*p))
  B <- c(tmp, -tmp) # diadic grid of univariate alternative values

  # find unique elements in the tail matrix, and establish a match
  unique_tail <- sort(unique(as.vector(tail)))
  tail_loc <- matrix(match(as.vector(tail), unique_tail), p)

  # update A and tail first
  A <- A + x_new
  tail <- tail + 1

  # expand A into a p x |B| matrix, and compute a p x |B| matrix R of tail loglik ratio stats
  A_expand <- matrix(A[cbind(rep(1:p, length(B)), as.vector(tail_loc))], p)
  R <- t(t(A_expand) * B - t(tail) * B^2 / 2)

  # reset some tail and CUSUM in A to zero
  tail <- tail * (R > 0)
  unique_tail_new <- sort(unique(as.vector(tail)))
  if (unique_tail_new[1] == 0){
    A <- A[, match(unique_tail_new[-1] - 1, unique_tail)]
    A <- cbind(0, A)
  } else {
    A <- A[, match(unique_tail_new - 1, unique_tail)]
  }

  # compute test stats S_diag, S_dense and S_sparse
  G <- sweep(A^2, 2, pmax(1, unique_tail_new), '/')  # normalisation
  G0 <- G * (G > a^2) # hard thresholded
  colsum_dense <- colSums(G)
  colsum_sparse <- colSums(G0)

  for (i in seq_along(unique_tail_new)) {
    colsum_dense[i] <- colsum_dense[i] - min(G[rowSums(tail == unique_tail_new[i]) > 0, i])
    colsum_sparse[i] <- colsum_sparse[i] - min(G0[rowSums(tail == unique_tail_new[i]) > 0, i])
  }

  S_diag <- max(R)
  S_dense <- max(colsum_dense)
  S_sparse <- max(colsum_sparse)

  stat <- setNames(c(S_diag, S_dense, S_sparse), c('diag', 'dense', 'sparse'))
  if (sparsity=='dense') stat <- stat[-3]
  if (sparsity=='sparse') stat <- stat[-2]

  return(list(stat=stat, A=A, tail=tail))
}

#' Processing a new data point for the 'Mei' class
#' @description This function implements the \code{\link{getData}} function to
#' perform the online changepoint detection for the 'Mei' class.
#' @param x_new a new data point
#' @param R vector of of tail CUSUMs that will be tracked and updated online
#' @param b a user specified lower bound on per-coordinate magnitude of the
#' vector of change to be detected.
#' @return a list of
#' \itemize{
#' \item stat: a vector of 2 test statistics for the 'Mei' class.
#' \item R: the updated R vector
#' }
#' @export
Mei_update <- function(x_new, R, b){
  R <- R + b * cbind(x_new, -x_new) - b^2/2  # update tail loglik ratio stats
  R[R<=0] <- 0 # reset when loglik ratio stat dips below 0

  stat <- setNames(c(max(R), max(colSums(R))), c('max', 'sum'))
  return(list(stat=stat,  R=R))
}


#' Processing a new data point for the 'Chan' class
#' @description This function implements the \code{\link{getData}} function to
#' perform the online changepoint detection for the 'Chan' class.
#' @param x_new a new data point
#' @param X_recent matrix of \code{w} most recent observations
#' @param CUSUM tail partial sums of different lengths to be tracked online
#' @param p0 sparsity parameter
#' @param w window parameter
#' @param lambda a tuning parameter for the 'Chan' method
#' @return a list of
#' \itemize{
#' \item stat: test statistic for the 'Chan' class.
#' \item X_recent: the updated X_recent matrix
#' \item CUSUM: the updated CUSUM matrix
#' }
#' @export
Chan_update <- function(x_new, X_recent, CUSUM, p0, w, lambda){
  CUSUM <- CUSUM + x_new - X_recent # update tail partial sums of length 1..w
  X_recent <- cbind(x_new, X_recent[, -w, drop=F])
  stat_plus <- max(colSums(log(1 - p0 + p0 * lambda * exp(t(t(pmax(CUSUM, 0))^2/(1:w))/4))))
  stat_minus <- max(colSums(log(1 - p0 + p0 * lambda * exp(t(t(pmax(-CUSUM, 0))^2/(1:w))/4))))
  return(list(stat=max(stat_plus, stat_minus),
              X_recent=X_recent,
              CUSUM=CUSUM))
}

#' Processing a new data point for the 'XS' class
#' @description This function implements the \code{\link{getData}} function to
#' perform the online changepoint detection for the 'XS' class.
#' @param x_new a new data point
#' @param X_recent matrix of \code{w} most recent observations
#' @param CUSUM tail partial sums of different lengths to be tracked online
#' @param p0 sparsity parameter
#' @param w window parameter
#' @return a list of
#' \itemize{
#' \item stat: test statistic for the 'XS' class.
#' \item X_recent: the updated X_recent matrix
#' \item CUSUM: the updated CUSUM matrix
#' }
#' @export
XS_update <- function(x_new, X_recent, CUSUM, p0, w){
  CUSUM <- CUSUM + x_new - X_recent # update tail partial sums of length 1..w
  X_recent <- cbind(x_new, X_recent[, -w, drop=F])
  stat_plus <- max(colSums(log(1 - p0 + p0 * exp(t(t(pmax(CUSUM, 0))^2/(1:w))/2))))
  stat_minus <- max(colSums(log(1 - p0 + p0 * exp(t(t(pmax(-CUSUM, 0))^2/(1:w))/2))))
  return(list(stat=max(stat_plus, stat_minus),
              X_recent=X_recent,
              CUSUM=CUSUM))
}


#' Compute Monte Carlo thresholds for the OCD method
#' @param dim Data dimension
#' @param patience Nominal patience of the procedure
#' @param beta Lower bound on l_2 norm of change
#' @param sparsity Sparsity parameter for the OCD method
#' @param MC_reps Number of Monte Carlo repetitions to use
#' @return A numeric vector of computed thresholds.
#' @export
MC_ocd <- function(dim, patience, beta, sparsity, MC_reps){
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c('diag','off_d','off_s')
  if (sparsity == 'sparse') peak_stat <- peak_stat[,-2]
  if (sparsity == 'dense') peak_stat <- peak_stat[,-3]

  # run MC_reps simulations for peak statistics of S_diag, S_{off,d} and S_{off,s}
  for (rep in 1:MC_reps){
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim))*2+4)

    for (i in 1:patience){
      x_new <- rnorm(dim)
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A; tail <- ret$tail
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}


#' Compute Monte Carlo thresholds for the Mei method
#' @param dim Data dimension
#' @param patience Nominal patience of the procedure
#' @param b lLwer bound on per-coordinate magnitude of change
#' @param MC_reps Number of Monte Carlo repetitions to use
#' @return A numeric vector of computed thresholds.
#' @export
MC_Mei <- function(dim, patience, b, MC_reps){
  peak_stat <- matrix(0, MC_reps, 2)
  colnames(peak_stat) <- c('max','sum')

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    R <- matrix(0, dim, 2)

    for (i in 1:patience){
      x_new <- rnorm(dim)
      ret <- Mei_update(x_new, R, b)
      R <- ret$R
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}

#' Compute Monte Carlo thresholds for the XS method
#' @param dim Data dimension
#' @param patience Nominal patience of the procedure
#' @param p0 Assumed fraction of nonzero coordinates of change.
#' @param w Window size
#' @param MC_reps number of Monte Carlo repetitions to use
#' @return A numeric vector of computed thresholds.
#' @export
MC_XS <- function(dim, patience, p0, w, MC_reps){
  peak_stat <- rep(-Inf, MC_reps)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:patience){
      x_new <- rnorm(dim)
      ret <- XS_update(x_new, X_recent, CUSUM, p0, w)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), exp(-1))
  return(th)
}

#' Compute Monte Carlo thresholds for the Chan method
#' @param dim Data dimension
#' @param patience Nominal patience of the procedure
#' @param p0 Assumed fraction of nonzero coordinates of change.
#' @param w Window size
#' @param lambda Tuning parameter for Chan (2017) mmethod
#' @param MC_reps number of Monte Carlo repetitions to use
#' @return A numeric vector of computed thresholds.
#' @export
MC_Chan <- function(dim, patience, p0, w, lambda, MC_reps){
  peak_stat <- rep(-Inf, MC_reps)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:patience){
      x_new <- rnorm(dim)
      ret <- Chan_update(x_new, X_recent, CUSUM, p0, w, lambda)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), exp(-1))
  return(th)
}
