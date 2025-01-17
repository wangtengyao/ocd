% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocd.R
\name{ChangepointDetector}
\alias{ChangepointDetector}
\title{Constructor for the ChangepointDetector S3 class}
\usage{
ChangepointDetector(dim, method = c("ocd", "Mei", "XS", "Chan"), thresh,
  patience = 5000, MC_reps = 100, beta = 1, sparsity = "auto",
  b = beta/sqrt(dim), p0 = 1/sqrt(dim), w = 200, lambda = sqrt(8) -
  2)
}
\arguments{
\item{dim}{Data dimension, all new data must be of this dimension}

\item{method}{Four methods are implemented: \code{ocd}, \code{Mei}, \code{XS}
and \code{Chan}. They correspond to the methods proposed in Chen, Wang and
Samworth (2020), Mei (2010), Xie and Siegmund (2013) and Chan (2017). The
constructed detector will be of 'OCD', 'Mei', 'XS' and 'Chan' subclass
respectively.}

\item{thresh}{A numeric vector or the character string 'MC'. If 'MC' is
specified then the correct threshold will be computed by Monte Carlo
simulation (the \code{patience} argument should be specified for this).
Otherwise, for method \code{ocd}, a vector of length 3 (corresponding
to the diagonal statistic, off-diagonal dense statistic and off-diagonal
sparse statistic) should be specifiied; for method \code{Mei}, a vector of
length two (corresponding to the max and sum statistics) should be specified;
for methods \code{XS} and \code{Chan}, a single positive real number should
be specified;}

\item{patience}{Required patience (average run length without change) of the
online changepoint procedure. This is optional if the thresholds for detection
are manually specified, but is required if Monte Carlo thresholds are used.}

\item{MC_reps}{Number of Monte Carlo repetitions to use to estimate the
thresholds. Only used when \code{thresh='MC'}.}

\item{beta}{lower bound on the l_2 norm of the vector of mean change to be
detected. This argument is used by the \code{ocd} method.}

\item{sparsity}{Parameter used by the \code{ocd}. If \code{sparsity='sparse'},
then only the diagonal and off-diagonal sparse statistics are used.
If \code{sparsity='dense'}, then only the diagonal and off-diagonal sparse
statistics are used. If \code{sparsity='auto'}, all three statistics are used
to detect both sparse and dense change adaptively.}

\item{b}{Lower bound on the per-coordinate magnitude of mean change be
detected. This argument is used by the 'Mei' method. If \code{b} is
unspecified but \code{beta} is specified, the default \code{b = beta/sqrt(dim)}
will be used.}

\item{p0}{A real number between 0 and 1. Sparsity parameter used by \code{XS}
and \code{Chan} methods. It is the assumed fraction of nonzero coordinates of
change. Default to \code{1/sqrt(dim)}.}

\item{w}{Window size parameter used by \code{XS} and \code{Chan} methods.
Number of most recent data points to keep track in memory. Default is 200.}

\item{lambda}{A tuning parameter used by the \code{Chan} method. Default is
\code{sqrt(8)-2}.}
}
\value{
An object of S3 class 'ChangepointDetector'. Depending on the
\code{method} argument specified, the object also belongs to a subclass
'OCD', 'Mei', 'XS' or 'Chan' corresponding to \code{method='ocd'}. It
contains the following attributes:
\itemize{
\item class - S3 class and subclass
\item data_dim - data dimension
\item method - method used for changepoint detection
\item param - a list of parameters used in the specific method: \code{beta}
and \code{sparsity} for method \code{ocd}; \code{b} for method \code{Mei};
\code{p0} and \code{w} for method \code{XS}; \code{p0}, \code{w} and
\code{lambda} for method \code{Chan}.
\item threshold - a named vector of thresholds used for detection (see the
\code{thresh} argument)
\item n_obs - number of observations, initialised to 0
\item baseline_mean - vector of pre-change mean, initialised to a vector of 0,
can be estimated by setting the changepoint detector into baseline mean and
standard deviation estimating status, see \code{\link{setStatus}}, or set
directly using \code{\link{setBaselineMean}}.
\item baseline_sd - vector of standard deviation, initialised to a vector of 1,
can be estimated by setting the changepoint detector into baseline mean and
standard deviation estimating status, see \code{\link{setStatus}}, or set
directly using \code{\link{setBaselineSD}}.
\item tracked -  a list of information tracked online by the changepoint
detector: matrices \code{A}
and \code{tail} for method \code{ocd}; vector \code{R} for method \code{Mei};
matrices \code{X_recent} and \code{CUSUM} for methods \code{XS} and \code{Chan}.
\item statistics - a named vector of test statistics for changepoint
detection: statistics with names \code{diag}, \code{off_d} and \code{off_s}
for method \code{ocd} (note if \code{sparsity} is \code{'dense'} or
\code{'sparse'}, then only (S^{diag}, S^{off,d})
and (S^{diag}, S^{off,s}) are included in \code{stat} respectively.);
statistics with names \code{max} and \code{sum} for
method \code{Mei}; a single numeric value for  methods \code{XS} and \code{Chan}.
\item status - one of the following: 'estimating' (the detector is estimating
the baseline mean and standard deviation with new data points), 'monitoring'
(the detector is detecting changes from the baseline mean from new data points)
and an integer recording the time of declaration of changepoint.
}
}
\description{
Constructor for the ChangepointDetector S3 class
}
\details{
This function is a wrapper. The \code{\link{new_OCD}},
\code{\link{new_Mei}}, \code{\link{new_XS}} and \code{\link{new_Chan}} carry
out the actual constructor implementation.
}
\examples{
detector_ocd <- ChangepointDetector(dim=100, method='ocd',
                                    thresh=c(11.6, 179.5, 54.9), beta=1)
detector_Mei <- ChangepointDetector(dim=100, method='Mei',
                                    thresh=c(8.6, 125.1), b=0.1)
detector_XS <- ChangepointDetector(dim=100, method='XS', thresh=55.1)
detector_Chan <- ChangepointDetector(dim=100, method='Chan', thresh=8.7)
}
\references{
\itemize{
\item Chen, Y., Wang, T. and Samworth, R. J. (2020) High-dimensional
multiscale online changepoint detection \emph{Preprint}. arxiv:2003.03668.
\item Mei, Y. (2010) Efficient scalable schemes for monitoring a large number
of data streams. \emph{Biometrika}, \strong{97}, 419--433.
\item Xie, Y. and Siegmund, D. (2013) Sequential multi-sensor change-point
detection.  \emph{Ann. Statist.}, \strong{41}, 670--692.
\item Chan, H. P. (2017) Optimal sequential detection in multi-stream data.
\emph{Ann. Statist.}, \strong{45}, 2736--2763.
}
}
\seealso{
accessor functions such as \code{\link{data_dim}}, the main function
for processing a new data point \code{\link{getData}}, other methods for the
ChangepointDetector class including \code{\link{reset}},
\code{\link{setBaselineMean}}, \code{\link{setBaselineSD}},
\code{\link{setStatus}}, \code{\link{normalisedStatistics}} and
\code{\link{checkChange}}.
}
