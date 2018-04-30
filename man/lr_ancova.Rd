\name{lr_ancova}
\alias{lr_ancova}

\title{Latent Regression for Group Effects with Latent Heteroskedastic
  Error Variance.}

\description{TBD}

\usage{
lr_ancova(outcome_model, Y, W, Z, G, varfuncs, plotfile="lr_ancova.pdf",
seed=12345, modelfileonly=FALSE, scalemat=NULL,blockprior=TRUE,...)
}

\arguments{
  
  \item{outcome_model}{TBD}

  \item{Y}{TBD}

  \item{W}{TBD}

  \item{Z}{TBD}

  \item{G}{TBD}

  \item{varfuncs}{TBD}

  \item{plotfile}{TBD}

  \item{seed}{TBD}

  \item{modelfileonly}{TBD}

  \item{scalemat}{TBD}

  \item{blockprior}{TBD}
  
  \item{...}{TBD}
}

\details{
  TBD
  Note that this function requires JAGS to be installed on your system.
}

\value{
  TBD
}

\references{
  Battauz, M. and Bellio, R. (2011). \dQuote{Structural modeling of
  measurement error in generalized linear models with Rasch measures
  as covariates,} \emph{Psychometrika} 76(1):40-56.
    
  Lockwood J.R. and McCaffrey D.F. (2014). \dQuote{Correcting for test
    score measurement error in ANCOVA models for estimating treatment
    effects,} \emph{Journal of Educational and Behavioral Statistics}
  39(1):22-52.

  Plummer, M. (2003). \dQuote{JAGS: A program for analysis of Bayesian
  graphical models using Gibbs sampling.} Proceedings of the 3rd
  International Workshop on Distributed Statistical Computing (DSC
  2003), Vienna, Austria.
}

\author{
  J.R. Lockwood \email{jrlockwood@ets.org}
}

\examples{
cat("TBD")
}

\keyword{models}
