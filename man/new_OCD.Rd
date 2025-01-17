% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocd.R
\name{new_OCD}
\alias{new_OCD}
\title{constructor of subclass 'OCD' in class 'ChangepointDetector'}
\usage{
new_OCD(dim, thresh, beta, sparsity)
}
\arguments{
\item{dim}{Data dimension, all new data must be of this dimension}

\item{thresh}{A numeric vector of length 3 (corresponding
to the diagonal statistic, off-diagonal dense statistic and off-diagonal
sparse statistic) should be specifiied.}

\item{beta}{Lower bound on the l_2 norm of the vector of mean change to be
detected.}

\item{sparsity}{If \code{sparsity='sparse'},
then only the diagonal and off-diagonal sparse statistics are used.
If \code{sparsity='dense'}, then only the diagonal and off-diagonal sparse
statistics are used. If \code{sparsity='auto'}, all three statistics are used
to detect both sparse and dense change adaptively.}
}
\value{
An object of S3 subclass 'OCD' in class 'ChangepointDetector'.
}
\description{
constructor of subclass 'OCD' in class 'ChangepointDetector'
}
\details{
It is preferred to use \code{\link{ChangepointDetector}} for
construction.
}
\examples{
detector <- new_OCD(dim=100, thresh=c(11.6, 179.5, 54.9), beta=1, sparsity='auto')
}
\references{
Chen, Y., Wang, T. and Samworth, R. J. (2020) High-dimensional multiscale
online changepoint detection \emph{Preprint}. arxiv:2003.03668.
}
