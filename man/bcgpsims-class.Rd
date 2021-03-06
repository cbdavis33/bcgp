% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{bcgpsims-class}
\alias{bcgpsims-class}
\title{An S4 class to represent a simulated model}
\description{
This class contains training and test data simulated directly from the BCGP
model along with some other information about the type of model the data was
simulated from (composite/non-composite, stationary/non-stationary,
deterministic/noisy). It will most commonly be useful after a call to
\code{simulate_from_model()} or \code{bcgpsims()}. The \code{bcgpsims} model
can then be plotted or sent to the function \code{bcgp_sampling} to fit the
model.
}
\section{Slots}{

\describe{
\item{\code{training}}{A list with elements \code{x}, an \code{n x d} matrix
containing the simulated training data locations, and \code{y}, a vector of
length \code{n} containing the simulated training data values. If the
\code{composite} argument was \code{TRUE} and the \code{decomposition}
argument was \code{TRUE} in the call to \code{simulate_from_model},
\code{training} will also have elements \code{yG},
\code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
\strong{L}ocal, and \strong{E}rror process values at the training data
locations.}

\item{\code{test}}{A list with elements \code{x}, an \code{nTest x d} matrix
containing the simulated test data locations, and \code{y}, a vector of
length \code{nTest} containing the simulated test data values. If
\code{composite} was \code{TRUE} and \code{decomposition} was \code{TRUE} in
the call to \code{simulate_from_model}, \code{test} will also have elements
\code{yG}, \code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
\strong{L}ocal, and \strong{E}rror process values at the test data locations.
Note: \code{yE} in \code{test} will \emph{always} be a zero vector. It is
only returned for completeness.}

\item{\code{parameters}}{A list containing the parameters used to simulate from the
model.}

\item{\code{stationary}}{A logical indicating whether the model is stationary or not.}

\item{\code{composite}}{A logical indicating whether the model is composite or not.}

\item{\code{seed}}{The seed used by the random number generator to generate the data.}
}}

\examples{
simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
}
\seealso{
\code{\link{simulate_from_model}} \code{\link{bcgpsims}}
}
