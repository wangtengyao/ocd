##### constructors #####
#' Constructor for the ChangepointDetector S3 class
#' @param dim Data dimension, all new data must be of this dimension
#' @param method Four methods are implemented: \code{ocd}, \code{Mei}, \code{XS}
#' and \code{Chan}. They correspond to the methods proposed in Chen, Wang and
#' Samworth (2020), Mei (2010), Xie and Siegmund (2013) and Chan (2017). The
#' constructed detector will be of 'OCD', 'Mei', 'XS' and 'Chan' subclass
#' respectively.
#' @param thresh A numeric vector or the character string 'MC'. If 'MC' is
#' specified then the correct threshold will be computed by Monte Carlo
#' simulation (the \code{patience} argument should be specified for this).
#' Otherwise, for method \code{ocd}, a vector of length 3 (corresponding
#' to the diagonal statistic, off-diagonal dense statistic and off-diagonal
#' sparse statistic) should be specifiied; for method \code{Mei}, a vector of
#' length two (corresponding to the max and sum statistics) should be specified;
#' for methods \code{XS} and \code{Chan}, a single positive real number should
#' be specified;
#' @param patience Required patience (average run length without change) of the
#' online changepoint procedure. This is optional if the thresholds for detection
#' are manually specified, but is required if Monte Carlo thresholds are used.
#' @param MC_reps Number of Monte Carlo repetitions to use to estimate the
#' thresholds. Only used when \code{thresh='MC'}.
#' @param beta lower bound on the l_2 norm of the vector of mean change to be
#' detected. This argument is used by the \code{ocd} method.
#' @param sparsity Parameter used by the \code{ocd}. If \code{sparsity='sparse'},
#' then only the diagonal and off-diagonal sparse statistics are used.
#' If \code{sparsity='dense'}, then only the diagonal and off-diagonal sparse
#' statistics are used. If \code{sparsity='auto'}, all three statistics are used
#' to detect both sparse and dense change adaptively.
#' @param b Lower bound on the per-coordinate magnitude of mean change be
#' detected. This argument is used by the 'Mei' method. If \code{b} is
#' unspecified but \code{beta} is specified, the default \code{b = beta/sqrt(dim)}
#' will be used.
#' @param p0 A real number between 0 and 1. Sparsity parameter used by \code{XS}
#' and \code{Chan} methods. It is the assumed fraction of nonzero coordinates of
#' change. Default to \code{1/sqrt(dim)}.
#' @param w Window size parameter used by \code{XS} and \code{Chan} methods.
#' Number of most recent data points to keep track in memory. Default is 200.
#' @param lambda A tuning parameter used by the \code{Chan} method. Default is
#' \code{sqrt(8)-2}.
#' @return An object of S3 class 'ChangepointDetector'. Depending on the
#' \code{method} argument specified, the object also belongs to a subclass
#' 'OCD', 'Mei', 'XS' or 'Chan' corresponding to \code{method='ocd'}. It
#' contains the following attributes:
#' \itemize{
#' \item class - S3 class and subclass
#' \item data_dim - data dimension
#' \item method - method used for changepoint detection
#' \item param - a list of parameters used in the specific method: \code{beta}
#' and \code{sparsity} for method \code{ocd}; \code{b} for method \code{Mei};
#' \code{p0} and \code{w} for method \code{XS}; \code{p0}, \code{w} and
#' \code{lambda} for method \code{Chan}.
#' \item threshold - a named vector of thresholds used for detection (see the
#' \code{thresh} argument)
#' \item n_obs - number of observations, initialised to 0
#' \item baseline_mean - vector of pre-change mean, initialised to a vector of 0,
#' can be estimated by setting the changepoint detector into baseline mean and 
#' standard deviation estimating status, see \code{\link{setStatus}}, or set 
#' directly using \code{\link{setBaselineMean}}.
#' \item baseline_sd - vector of standard deviation, initialised to a vector of 1, 
#' can be estimated by setting the changepoint detector into baseline mean and 
#' standard deviation estimating status, see \code{\link{setStatus}}, or set 
#' directly using \code{\link{setBaselineSD}}.
#' \item tracked -  a list of information tracked online by the changepoint
#' detector: matrices \code{A}
#' and \code{tail} for method \code{ocd}; vector \code{R} for method \code{Mei};
#' matrices \code{X_recent} and \code{CUSUM} for methods \code{XS} and \code{Chan}.
#' \item statistics - a named vector of test statistics for changepoint
#' detection: statistics with names \code{diag}, \code{off_d} and \code{off_s}
#' for method \code{ocd} (note if \code{sparsity} is \code{'dense'} or
#' \code{'sparse'}, then only (S^{diag}, S^{off,d})
#' and (S^{diag}, S^{off,s}) are included in \code{stat} respectively.);
#' statistics with names \code{max} and \code{sum} for
#' method \code{Mei}; a single numeric value for  methods \code{XS} and \code{Chan}.
#' \item status - one of the following: 'estimating' (the detector is estimating
#' the baseline mean and standard deviation with new data points), 'monitoring' 
#' (the detector is detecting changes from the baseline mean from new data points) 
#' and an integer recording the time of declaration of changepoint.
#' }
#' @details This function is a wrapper. The \code{\link{new_OCD}},
#' \code{\link{new_Mei}}, \code{\link{new_XS}} and \code{\link{new_Chan}} carry
#' out the actual constructor implementation.
#' @seealso accessor functions such as \code{\link{data_dim}}, the main function
#' for processing a new data point \code{\link{getData}}, other methods for the
#' ChangepointDetector class including \code{\link{reset}},
#' \code{\link{setBaselineMean}}, \code{\link{setBaselineSD}}, 
#' \code{\link{setStatus}}, \code{\link{normalisedStatistics}} and 
#' \code{\link{checkChange}}.
#' @examples
#' detector_ocd <- ChangepointDetector(dim=100, method='ocd',
#'                                     thresh=c(11.6, 179.5, 54.9), beta=1)
#' detector_Mei <- ChangepointDetector(dim=100, method='Mei',
#'                                     thresh=c(8.6, 125.1), b=0.1)
#' detector_XS <- ChangepointDetector(dim=100, method='XS', thresh=55.1)
#' detector_Chan <- ChangepointDetector(dim=100, method='Chan', thresh=8.7)
#' @references
#' \itemize{
#' \item Chen, Y., Wang, T. and Samworth, R. J. (2020) High-dimensional
#' multiscale online changepoint detection \emph{Preprint}. arxiv:2003.03668.
#' \item Mei, Y. (2010) Efficient scalable schemes for monitoring a large number
#' of data streams. \emph{Biometrika}, \strong{97}, 419--433.
#' \item Xie, Y. and Siegmund, D. (2013) Sequential multi-sensor change-point
#' detection.  \emph{Ann. Statist.}, \strong{41}, 670--692.
#' \item Chan, H. P. (2017) Optimal sequential detection in multi-stream data.
#' \emph{Ann. Statist.}, \strong{45}, 2736--2763.
#' }
#' @export
ChangepointDetector <- function(dim, method=c('ocd', 'Mei', 'XS', 'Chan'),
                                thresh, patience=5000, MC_reps=100,
                                beta=1, sparsity='auto', b=beta/sqrt(dim),
                                p0=1/sqrt(dim), w=200, lambda=sqrt(8)-2){
  if (identical(thresh, 'MC')){
    thresh <- switch(method,
                     ocd = MC_ocd(dim, patience, beta, sparsity, MC_reps),
                     Mei = MC_Mei(dim, patience, b, MC_reps),
                     XS = MC_XS(dim, patience, p0, w, MC_reps),
                     Chan = MC_Chan(dim, patience, p0, w, lambda, MC_reps))
  }
  detector <- switch(method,
         ocd = new_OCD(dim, thresh, beta, sparsity),
         Mei = new_Mei(dim, thresh, b),
         XS = new_XS(dim, thresh, p0, w),
         Chan = new_Chan(dim, thresh, p0, w, lambda))
  return(detector)
}

#' constructor of subclass 'OCD' in class 'ChangepointDetector'
#' @param dim Data dimension, all new data must be of this dimension
#' @param thresh A numeric vector of length 3 (corresponding
#' to the diagonal statistic, off-diagonal dense statistic and off-diagonal
#' sparse statistic) should be specifiied.
#' @param beta Lower bound on the l_2 norm of the vector of mean change to be
#' detected.
#' @param sparsity If \code{sparsity='sparse'},
#' then only the diagonal and off-diagonal sparse statistics are used.
#' If \code{sparsity='dense'}, then only the diagonal and off-diagonal sparse
#' statistics are used. If \code{sparsity='auto'}, all three statistics are used
#' to detect both sparse and dense change adaptively.
#' @return An object of S3 subclass 'OCD' in class 'ChangepointDetector'.
#' @details It is preferred to use \code{\link{ChangepointDetector}} for
#' construction.
#' @examples
#' detector <- new_OCD(dim=100, thresh=c(11.6, 179.5, 54.9), beta=1, sparsity='auto')
#' @references
#' Chen, Y., Wang, T. and Samworth, R. J. (2020) High-dimensional multiscale
#' online changepoint detection \emph{Preprint}. arxiv:2003.03668.
#' @export
new_OCD <- function(dim, thresh, beta, sparsity){
  L <- floor(log2(dim))*2+4
  A <- matrix(0, dim, 1)
  tail <- matrix(0, dim, L)
  stats <- setNames(c(0,0,0), c('diag','off_d','off_s'))
  if (sparsity=='sparse') stats <- stats[-2]
  if (sparsity=='dense') stats <- stats[-3]

  detector <- structure(list(),
                        class = c('OCD', 'ChangepointDetector'),
                        data_dim = dim,                    # immutable attributes
                        method = 'ocd',
                        param = list(beta=beta, sparsity=sparsity),
                        thresholds = thresh,
                        n_obs = 0,                          # mutable quantities
                        baseline_mean = rep(0,dim),
                        baseline_sd = rep(1,dim),
                        tracked = list(A=A, tail=tail),
                        statistics = stats,
                        status = 'monitoring')
  return(detector)
}

#' constructor of subclass 'Mei' in class 'ChangepointDetector'
#' @param dim Data dimension, all new data must be of this dimension
#' @param thresh Detection threshold. A numeric vector of
#' length two (corresponding to the max and sum statistics) should be specified.
#' @param b Lower bound on the per-coordinate magnitude of mean change be
#' detected.
#' @return An object of S3 subclass 'Mei' in class 'ChangepointDetector'.
#' @details It is preferred to use \code{\link{ChangepointDetector}} for
#' construction.
#' @examples
#' detector <- new_Mei(dim=100, thresh=c(8.6, 125.1), b=0.1)
#' @references
#' Mei, Y. (2010) Efficient scalable schemes for monitoring a large number
#' of data streams. \emph{Biometrika}, \strong{97}, 419--433.
#' @export
new_Mei <- function(dim, thresh, b){
  R <- matrix(0, dim, 2)
  stats <- setNames(c(0,0), c('max','sum'))

  detector <- structure(list(),
                        class = c('Mei', 'ChangepointDetector'),
                        data_dim = dim,
                        method = 'Mei',
                        param = list(b=b),
                        thresholds = thresh,
                        n_obs = 0,
                        baseline_mean = rep(0, dim),
                        baseline_sd = rep(1,dim),
                        tracked = list(R=R),
                        statistics = stats,
                        status = 'monitoring')
  return(detector)
}

#' constructor of subclass 'XS' in class 'ChangepointDetector'
#' @param dim Data dimension, all new data must be of this dimension
#' @param thresh Detection threshold. A positive real number.
#' @param p0 A sparsity parameter between 0 and 1. It is the assumed fraction of
#' nonzero coordinates of change. Default to \code{1/sqrt(dim)}.
#' @param w Window size parameter.
#' Number of most recent data points to keep track in memory. Default is 200.
#' @return An object of S3 subclass 'XS' in class 'ChangepointDetector'.
#' @details It is preferred to use \code{\link{ChangepointDetector}} for
#' construction.
#' @examples
#' detector <- new_XS(dim=100, thresh=55.1, p0=0.1, w=200)
#' @references
#' Xie, Y. and Siegmund, D. (2013) Sequential multi-sensor change-point
#' detection.  \emph{Ann. Statist.}, \strong{41}, 670--692.
#' @export
new_XS <- function(dim, thresh, p0, w){
  X_recent <- matrix(0, dim, w)
  CUSUM <- matrix(0, dim, w)
  stats <- setNames(0, "")

  detector <- structure(list(),
                        class = c('XS', 'ChangepointDetector'),
                        data_dim = dim,
                        method = 'XS',
                        param = list(p0=p0, w=w),
                        thresholds = thresh,
                        n_obs = 0,
                        baseline_mean = rep(0, dim),
                        baseline_sd = rep(1,dim),
                        tracked = list(X_recent=X_recent, CUSUM=CUSUM),
                        statistics = stats,
                        status = 'monitoring')
  return(detector)
}

#' construtor for subclass 'Chan' in class 'ChangepointDetector'
#' @param dim Data dimension, all new data must be of this dimension
#' @param thresh Detection threshold. A positive real number.
#' @param p0 A sparsity parameter between 0 and 1. It is the assumed fraction of
#' nonzero coordinates of change. Default to \code{1/sqrt(dim)}.
#' @param w Window size parameter.
#' Number of most recent data points to keep track in memory. Default is 200.
#' @param lambda A tuning parameter used by the Chan (2017) method.  Default is
#' \code{sqrt(8)-2}.
#' @return An object of S3 subclass 'Chan' in class 'ChangepointDetector'.
#' @details It is preferred to use \code{\link{ChangepointDetector}} for
#' construction.
#' @examples
#' detector <- new_Chan(dim=100, thresh=8.7, p0=0.1, w=200, lambda=sqrt(8)-2)
#' @references
#' Chan, H. P. (2017) Optimal sequential detection in multi-stream data.
#' \emph{Ann. Statist.}, \strong{45}, 2736--2763.
#' @export
new_Chan <- function(dim, thresh, p0, w, lambda){
  X_recent <- matrix(0, dim, w)
  CUSUM <- matrix(0, dim, w)
  stats <- setNames(0, "")

  detector <- structure(list(),
                        class = c('Chan', 'ChangepointDetector'),
                        data_dim = dim,
                        method = 'Chan',
                        param = list(p0=p0, w=w, lambda=lambda),
                        thresholds = thresh,
                        n_obs = 0,
                        baseline_mean = rep(0, dim),
                        baseline_sd = rep(1,dim),
                        tracked = list(X_recent=X_recent, CUSUM=CUSUM),
                        statistics = stats,
                        status = 'monitoring')
  return(detector)
}

##### methods for the class ChangepointDetector and its subclasses #####

#' Accessor functions to attributes of class ChangepointDetector
#' @param detector Object of S3 class 'ChangepointDetector'
#' @details Obtain various attributes of the class 'ChangepointDetector'
#' @seealso \code{\link{ChangepointDetector}}
#' @name accessor
NULL

#' @rdname accessor
#' @export
data_dim <- function(detector) attr(detector, 'data_dim')

#' @rdname accessor
#' @export
ocdMethod <- function(detector) attr(detector, 'method')

#' @rdname accessor
#' @export
n_obs <- function(detector) attr(detector, 'n_obs')

#' @rdname accessor
#' @export
patience <- function(detector) attr(detector, 'patience')

#' @rdname accessor
#' @export
param <- function(detector) attr(detector, 'param')

#' @rdname accessor
#' @export
thresholds <- function(detector) attr(detector, 'thresholds')

#' @rdname accessor
#' @export
baselineMean <- function(detector) attr(detector, 'baseline_mean')

#' @rdname accessor
#' @export
baselineSD <- function(detector) attr(detector, 'baseline_sd')

#' @rdname accessor
#' @export
tracked <- function(detector) attr(detector, 'tracked')

#' @rdname accessor
#' @export
statistics <- function(detector) attr(detector, 'statistics')

#' @rdname accessor
#' @export
status <- function(detector) attr(detector, 'status')


#' Reset changepoint detector to initial state
#' @param detector Object of class 'Changepoint Detector'
#' @return Updated object \code{detector}
#' @export
reset <- function(detector) UseMethod('reset')


#' @describeIn reset Reset object of subclass 'OCD'
#' @export
reset.OCD <- function(detector){
  p <- data_dim(detector)
  L <- floor(log2(p))*2+4
  A <- matrix(0, p, 1)
  tail <- matrix(0, p, L)
  stats <- setNames(c(0,0,0), c('diag','off_d','off_s'))
  if (param(detector)$sparsity=='sparse') stats <- stats[-2]
  if (param(detector)$sparsity=='dense') stats <- stats[-3]

  attr(detector, 'n_obs') <- 0
  attr(detector, 'baseline_mean') <- rep(0,p)
  attr(detector, 'baseline_sd') <- rep(1,p)
  attr(detector, 'tracked') <- list(A=A, tail=tail)
  attr(detector, 'statistics') <- stats
  attr(detector, 'status') <- 'monitoring'
  return(detector)
}

#' @describeIn reset Reset object of subclass 'Mei'
#' @export
reset.Mei <- function(detector){
  p <- data_dim(detector)
  R <- matrix(0, p, 2)
  stats <- setNames(c(0,0), c('max','sum'))

  attr(detector, 'n_obs') <- 0
  attr(detector, 'baseline_mean') <- rep(0,p)
  attr(detector, 'baseline_sd') <- rep(1,p)
  attr(detector, 'tracked') <- list(R=R)
  attr(detector, 'statistics') <- stats
  attr(detector, 'status') <- 'monitoring'
  return(detector)
}

#' @describeIn reset Reset object of subclass 'XS'
#' @export
reset.XS <- function(detector){
  p <- data_dim(detector)
  tmp <- param(detector); w <- tmp$w
  X_recent <- matrix(0, p, w)
  CUSUM <- matrix(0, p, w)
  stats <- 0

  attr(detector, 'n_obs') <- 0
  attr(detector, 'baseline_mean') <- rep(0,p)
  attr(detector, 'baseline_sd') <- rep(1,p)
  attr(detector, 'tracked') <- list(X_recent=X_recent, CUSUM=CUSUM)
  attr(detector, 'statistics') <- stats
  attr(detector, 'status') <- 'monitoring'
  return(detector)
}

#' @describeIn reset Reset object of subclass 'Chan'
#' @export
reset.Chan <- function(detector){
  p <- data_dim(detector)
  tmp <- param(detector); w <- tmp$w
  X_recent <- matrix(0, p, w)
  CUSUM <- matrix(0, p, w)
  stats <- 0

  attr(detector, 'n_obs') <- 0
  attr(detector, 'baseline_mean') <- rep(0,p)
  attr(detector, 'baseline_sd') <- rep(1,p)
  attr(detector, 'tracked') <- list(X_recent=X_recent, CUSUM=CUSUM)
  attr(detector, 'statistics') <- stats
  attr(detector, 'status') <- 'monitoring'
  return(detector)
}


#' Set baseline mean
#' @param detector Object of class 'Changepoint Detector'
#' @param mean vector of pre-change mean, must be of the same dimension as
#' specified in the data_dim attribute of \code{detector}.
#' @return Updated object \code{detector}
#' @export
setBaselineMean <- function(detector, mean){
  attr(detector, 'baseline_mean') <- mean
  return(detector)
}

#' Set baseline standard deviation
#' @param detector Object of class 'Changepoint Detector'
#' @param sd vector of standard deviation, must be of the same dimension as
#' specified in the data_dim attribute of \code{detector}.
#' @return Updated object \code{detector}
#' @export
setBaselineSD <- function(detector, sd){
  attr(detector, 'baseline_sd') <- sd
  return(detector)
}

#' Set changepoint detector status
#' @param detector Object of class 'Changepoint Detector'
#' @param new_status 'estimating' or 'monitoring'
#' @return Updated object \code{detector}
#' @details If the status is set to 'estimating', new observations are used to
#' update current estimate of pre-change mean and standard deviation. If the 
#' status is set to 'monitoring', new observations are used to check if mean 
#' change has occurred.
#' @export
setStatus <- function(detector, new_status){
  attr(detector, 'status') <- new_status
  if (new_status=='estimating'){
    detector <- setBaselineMean(detector, rep(0, data_dim(detector)))
    detector <- setBaselineSD(detector, rep(1, data_dim(detector)))
  } else if (new_status=='monitoring'){
    attr(detector, 'n_obs') <- 0   # reset observation counter for detection
  }
  return(detector)
}

#' Compute maximum ratio between detection statistic and its threshold
#' @param detector Object of class 'Changepoint Detector'
#' @return maximum of the ratio between the current test statistics and their
#' respective thresholds.
#' @export
normalisedStatistics <- function(detector) {
  max(statistics(detector) / thresholds(detector))
}

#' Check if a mean change has occurred.
#' @param detector Object of class 'Changepoint Detector'
#' @return Updated object \code{detector}
#' @details The \code{\link{normalisedStatistics}} funcrtion is used to check if
#' any of the test statistics are above the threshold level. If this happens, the
#' status of the detector is changed to record the time of change and a message
#' is printed to the standard output declaring the change.
#' @seealso \code{\link{normalisedStatistics}}, \code{\link{setStatus}},
#' @export
checkChange <- function(detector){
  if (normalisedStatistics(detector) >= 1){
    n <- n_obs(detector)
    attr(detector, 'status') <- setNames(n, 'declared at')
    cat('Changepoint declared at time =', n, '\n')
  }
  return(detector)
}


#' Processing a new data point
#' @description This is the main function for the 'ChangepointDetector' class.
#' @param detector Object of class 'Changepoint Detector'
#' @param x_new A new data point. It must be of the same dimension as
#' specified in the data_dim attribute of \code{detector}.
#' @return Updated object \code{detector}
#' @details If the status of the \code{detector} object is 'estimating', the new
#' data point is used to update the current estimate of pre-change mean and 
#' standard deviation. If the status of the \code{detector} object is monitoring', 
#' the new data point is used to detect if a mean change has occurred.
#' @seealso \code{\link{setBaselineMean}} for updating the pre-change mean
#' estimate, \code{\link{setBaselineSD}} for updating the standard deviation 
#' estimate, \code{\link{checkChange}} for checking for change.
#' @export
getData <- function(detector, x_new) UseMethod('getData')


#' @describeIn getData Process a new data for subclass 'OCD'
#' @export
getData.OCD <- function(detector, x_new){
  attr(detector, 'n_obs') <- attr(detector, 'n_obs') + 1
  if (status(detector)=='estimating'){  # use the new data to update baseline mean and standard deviation estimates
    if (n_obs(detector) == 1){
      new_mean <- x_new
      new_sd <- 0
    } else {
      new_mean <- (x_new + baselineMean(detector) * (n_obs(detector) - 1)) / n_obs(detector)
      new_sd <- sqrt(((n_obs(detector) - 2) * baselineSD(detector)^2 + (n_obs(detector) - 1) * baselineMean(detector)^2 + 
                        x_new^2 - n_obs(detector) * new_mean^2) / (n_obs(detector) - 1))
    }
    detector <- setBaselineMean(detector, new_mean)
    detector <- setBaselineSD(detector, new_sd)
  } else { # use the new data to update tracked info and compute test stats
    tmp <- tracked(detector); A <- tmp$A; tail <- tmp$tail
    tmp <- param(detector); beta <- tmp$beta; sparsity <- tmp$sparsity
    ret <- ocd_update((x_new - baselineMean(detector))/baselineSD(detector), A, tail, beta, sparsity)
    attr(detector, 'statistics') <- ret$stat
    attr(detector, 'tracked') <- list(A=ret$A, tail=ret$tail)
    if (status(detector)=='monitoring') detector <- checkChange(detector)
  }
  return(detector)
}

#' @describeIn getData Process a new data for subclass 'Mei'
#' @export
getData.Mei <- function(detector, x_new){
  attr(detector, 'n_obs') <- attr(detector, 'n_obs') + 1
  if (status(detector)=='estimating'){  # use the new data to update baseline mean estimate
    if (n_obs(detector) == 1){
      new_mean <- x_new
      new_sd <- 0
    } else {
      new_mean <- (x_new + baselineMean(detector) * (n_obs(detector) - 1)) / n_obs(detector)
      new_sd <- sqrt(((n_obs(detector) - 2) * baselineSD(detector) + (n_obs(detector) - 1) * baselineMean(detector)^2 + 
                        x_new^2 - n_obs(detector) * new_mean^2) / (n_obs(detector) - 1))
    }
    detector <- setBaselineMean(detector, new_mean)
    detector <- setBaselineSD(detector, new_sd)
  } else { # use the new data to update tracked info and compute test stats
    tmp <- tracked(detector); R <- tmp$R
    tmp <- param(detector); b <- tmp$b
    ret <- Mei_update((x_new - baselineMean(detector))/baselineSD(detector), R, b)
    attr(detector, 'statistics') <- ret$stat
    attr(detector, 'tracked') <- list(R=ret$R)
    if (status(detector)=='monitoring') detector <- checkChange(detector)
  }
  return(detector)
}

#' @describeIn getData Process a new data for subclass 'XS'
#' @export
getData.XS <- function(detector, x_new){
  attr(detector, 'n_obs') <- attr(detector, 'n_obs') + 1
  if (status(detector)=='estimating'){  # use the new data to update baseline mean estimate
    if (n_obs(detector) == 1){
      new_mean <- x_new
      new_sd <- 0
    } else {
      new_mean <- (x_new + baselineMean(detector) * (n_obs(detector) - 1)) / n_obs(detector)
      new_sd <- sqrt(((n_obs(detector) - 2) * baselineSD(detector) + (n_obs(detector) - 1) * baselineMean(detector)^2 + 
                        x_new^2 - n_obs(detector) * new_mean^2) / (n_obs(detector) - 1))
    }
    detector <- setBaselineMean(detector, new_mean)
    detector <- setBaselineSD(detector, new_sd)
  } else { # use the new data to update tracked info and compute test stats
    tmp <- tracked(detector); X_recent <- tmp$X_recent; CUSUM <- tmp$CUSUM
    tmp <- param(detector); p0 <- tmp$p0; w <- tmp$w
    ret <- XS_update((x_new - baselineMean(detector))/baselineSD(detector), X_recent, CUSUM, p0, w)
    attr(detector, 'statistics') <- ret$stat
    attr(detector, 'tracked') <- list(X_recent=ret$X_recent, CUSUM=ret$CUSUM)
    if (status(detector)=='monitoring') detector <- checkChange(detector)
  }
  return(detector)
}

#' @describeIn getData Process a new data for subclass 'Chan'
#' @export
getData.Chan <- function(detector, x_new){
  attr(detector, 'n_obs') <- attr(detector, 'n_obs') + 1
  if (status(detector)=='estimating'){  # use the new data to update baseline mean estimate
    if (n_obs(detector) == 1){
      new_mean <- x_new
      new_sd <- rep(0,length(x_new))
    } else {
      new_mean <- (x_new + baselineMean(detector) * (n_obs(detector) - 1)) / n_obs(detector)
      new_sd <- sqrt(((n_obs(detector) - 2) * baselineSD(detector) + (n_obs(detector) - 1) * baselineMean(detector)^2 + 
                        x_new^2 - n_obs(detector) * new_mean^2) / (n_obs(detector) - 1))
    }
    detector <- setBaselineMean(detector, new_mean)
    detector <- setBaselineSD(detector, new_sd)
  } else { # use the new data to update tracked info and compute test stats
    tmp <- tracked(detector); X_recent <- tmp$X_recent; CUSUM <- tmp$CUSUM
    tmp <- param(detector); p0 <- tmp$p0; w <- tmp$w; lambda <- tmp$lambda
    ret <- Chan_update((x_new - baselineMean(detector))/baselineSD(detector), X_recent, CUSUM, p0, w, lambda)
    attr(detector, 'statistics') <- ret$stat
    attr(detector, 'tracked') <- list(X_recent=ret$X_recent, CUSUM=ret$CUSUM)
    if (status(detector)=='monitoring') detector <- checkChange(detector)
  }
  return(detector)
}


#' Printing methods for the 'ChangepointDetector' class
#' @param x object of the 'ChangepointDetector' class
#' @param ... other arguments used in \code{print}
#' @export
print.ChangepointDetector <- function(x, ...){
  detector <- x
  cat('Online changepoint detector using method:', ocdMethod(detector), '\n\n')
  cat('Time =', n_obs(detector), '\n\n')
  baseline_mean <- round(baselineMean(detector), 3)
  baseline_sd <- round(baselineSD(detector), 3)
  if (data_dim(detector) <= 10){
    cat('Baseline mean =', baseline_mean, '\n')
    cat('Standard deviation =', baseline_sd, '\n\n')
  } else {
    cat('Baseline mean =', head(baseline_mean,8), '...',
    tail(baseline_mean,2), '\n')
    cat('Standard deviation =', head(baseline_sd,8), '...', 
    tail(baseline_sd,2), '\n\n')
  }

  if (status(detector)!='estimating'){
    mx <- rbind(statistics(detector), thresholds(detector))
    mx <- round(mx, 3)
    row.names(mx) <- c('statistics', 'thresholds')
    print(mx)
  }

  if (is.numeric(status(detector))){
    cat('\nChangepoint declared at time =', status(detector), '\n\n')
  }
}
