\name{triple_goal}
\alias{triple_goal}

\title{Shen and Louis (1998) Triple Goal Estimators}

\description{ \code{triple_goal} implements the "Triple Goal" estimates
  of Shen and Louis (1998) for a vector of parameters given a sample
  from the posterior distribution of those parameters.  Also computes
  "constrained Bayes" estimators of Ghosh (1992).}

\usage{
triple_goal(s, stop.if.ties = FALSE, quantile.type = 7)
}

\arguments{
  \item{s}{
    A (n x K) matrix of "n" samples of "K" group parameters.
  }
  
  \item{stop.if.ties}{
    logical; if TRUE, function stops if any units have identical
    posterior mean ranks; otherwise breaks ties at random.
  }
  
  \item{quantile.type}{
    "type" argument to \code{\link{quantile}} function for different
    methods of computing quantiles.
  }
}

\details{
  In typical applications, the matrix \code{s} will be a sample of size
  "n" from the joint posterior distribution of a vector of
  "K" group-specific parameters.  Both the triple goal and constrained
  Bayes estimators are designed to mitigate problems arising from
  underdispersion of posterior means; see references.
}

\value{A dataframe with "K" rows with fields:
  \item{theta_pm}{Posterior mean estimates of group parameters.}

  \item{theta_psd}{Posterior standard deviation estimates of group
    parameters.}

  \item{theta_cb}{"Constrained Bayes" estimates of group parameters
    using formula in Shen and Louis (1998).}

  \item{theta_gr}{"Triple Goal" estimates of group parameters using
    algorithm defined in Shen and Louis (1998).}

  \item{rbar}{Posterior means of ranks of group parameters (1=lowest).}

  \item{rhat}{Integer ranks of group parameters (=rank(rbar)).}
}

\references{
  Wei Shen and Thomas A. Louis (1998).
  \dQuote{Triple-goal estimates in two-stage hierarchical models,}
  \emph{Journal of the Royal Statistical Society, Series B} 60(2):455-471.

  Malay Ghosh (1992).
  \dQuote{Constrained Bayes estimation with applications,}
  \emph{Journal of the American Statistical Association} 87(418):533-540.
}

\author{
  J.R. Lockwood \email{jrlockwood@ets.org}
}

\examples{
.K <- 50
.nsamp <- 500
.theta_true <- rnorm(.K)
.s <- matrix(.theta_true, ncol=.K, nrow=.nsamp, byrow=TRUE) + matrix(rnorm(.K*.nsamp, sd=0.4), ncol=.K, nrow=.nsamp)
.e <- triple_goal(.s)
str(.e)
head(.e)
}

\keyword{utilities}
