#' ocd: A package for high-dimensional multiscale online changepoint detection
#'
#' The ocd package provides the S3 class \code{\link{ChangepointDetector}} that
#' processes data sequentially using the \code{\link{getData}} function and
#' aims to detect change as soon as it occurs online subject to false alarm
#' rates.
#'
#' @seealso \code{\link{ChangepointDetector}} for detailed usage.
#' @references Chen, Y., Wang, T. and Samworth, R. J. (2020) High-dimensional
#' multiscale online changepoint detection \emph{Preprint}.
#'
#' @docType package
#' @name ocd
#' @examples
#' library(magrittr)
#' p <- 100
#' thresh <- setNames(c(11.62, 179.48, 54.87), c('diag', 'off_d', 'off_s'))
#' detector <- new_OCD(p=p, beta=1, thresh=thresh, sparsity='auto')
#' detector %<>% setStatus('estimating')
#' for (i in 1:10000){
#'   x_new <- rnorm(p, mean=old_mean)
#'   detector %<>% getData(x_new)
#' }
#' detector %>% print
#'
#' detector %<>% setStatus('monitoring')
#' for (i in 1:500){
#'   x_new <- rnorm(p, old_mean * (i < 250) + new_mean * (i>=250))
#'   detector %<>% getData(x_new)
#' }
#' detector %>% print
NULL
