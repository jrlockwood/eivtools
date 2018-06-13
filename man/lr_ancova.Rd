\name{lr_ancova}
\alias{lr_ancova}

\title{Latent Regression for Group Effects with Latent Heteroskedastic
  Error Variance
}

\description{Uses the \code{\link{jags}} function in R2jags to fit a
  latent-variable GLM with error-prone covariates that may have
  heteroskedastic normal measurement error with variance that is a
  function of the latent variable, such as commonly occurs with test
  scores computed using item-response-theory (IRT) models.
}

\usage{
lr_ancova(outcome_model, Y, W, Z, G, varfuncs, plotfile=NULL,
seed=12345, modelfileonly=FALSE, scalemat=NULL, blockprior=TRUE, ...)
}

\arguments{
  
  \item{outcome_model}{A character string indicating the outcome model.
    Valid values are currently 'normal', 'poisson',
    'bernoulli_probit', and 'bernoulli_logit'.}

  \item{Y}{A numeric vector of outcome values. Missing values are allowed.}

  \item{W}{A numeric matrix of error-prone covariates.  Missing values
    are allowed, though no column of \code{W} may be entirely missing.}

  \item{Z}{A numeric matrix of error-free covariates.  Missing values
    are not allowed.  First column must be a vector of 1s to serve as a
    model intercept because effects of groups in \code{G} are
    parameterized with a sum-to-zero constraint.  See Details for
    additional information.}

  \item{G}{A numeric or factor vector indicating group memberships of
    units.  Missing values not allowed.}

  \item{varfuncs}{A list with as many components as there are
    error-prone covariates, equal to the number of columns of \code{W}.
    For each \code{i} from \code{1} to \code{ncol(W)},
    \code{varfuncs[[i]]} is itself a list summarizing the known
    information about the measurement error in the variable
    \code{W[,i]}.  See Details.}

  \item{plotfile}{Character string providing full path to a PDF file
    that will store some diagnostic plots regarding the variance
    functions.  Default is NULL and will be assigned to a file in a
    temporary directory and the name of file will be returned.}

  \item{seed}{An integer that will be passed to \code{\link{set.seed}}
    so that Monte Carlo results can be reproduced.}

  \item{modelfileonly}{If TRUE, function will return a link to a file
    that contains the JAGS model code, but will not actually fit the
    model.  Default is FALSE.}

  \item{scalemat}{When there are multiple error-prone covariates, the
    specification of the Bayesian model as implemented in JAGS requires
    a scale matrix for a Wishart prior distribution for a precision
    matrix.  The default is NULL, in which case the function will set a
    value of \code{scalemat}; see Details.  If the user wishes to pass a
    \code{scalemat} it must be a positive-definite symmetric matrix of
    dimension \code{ncol(W)}.}

  \item{blockprior}{If TRUE (the default), specifies JAGS code to
    encourage updating regression model parameters as a block to improve
    MCMC mixing.}

  \item{ ... }{Additional arguments to \code{\link{jags}}.}
}

\details{

  \emph{Theory}

  The outcome \eqn{Y}{Y} is assumed to depend on \eqn{(X,Z,G)}{(X,Z,G)}
  where \eqn{X}{X} is a vector of latent variables, \eqn{Z}{Z} is a
  vector of observed, error-free variables, and \eqn{G}{G} is a grouping
  variable.  For example, one may be interested in the effects of some
  intervention where \eqn{G}{G} indicates groupings of units that
  received different treatments, and the variables \eqn{(X,Z)}{(X,Z)}
  are potential confounders.  This function addresses the case where
  \eqn{X}{X} is unobserved, and error-prone proxies \eqn{W}{W} are
  instead observed.  It is assumed that \eqn{W = X + U}{W = X + U} for
  mean-zero, normally-distributed measurement error \eqn{U}{U}, and that
  \eqn{Var(U)} may be a function \eqn{g(X)}{g(X)} of \eqn{X}{X}.  Such
  error structures commonly arise with the use of test scores computed
  using item-response-theory (IRT) models. Details on these issues and
  other model assumptions are provided in the references.  The model is
  a generalization of errors-in-variables linear regression.

  The model assumes that the outcome \eqn{Y}{Y} depends on
  \eqn{(X,Z,G)}{(X,Z,G)} through a linear function of these predictors,
  and parameters for this linear function are estimated.  The
  conditional distribution of \eqn{Y}{Y} given these predictors may be
  normal, Poisson, or Bernoulli.  If normal, the model also estimates a
  residual variance for \eqn{Y}{Y} given the covariates.  If Poisson, a
  log link is assumed.  If Bernoulli, either a logit or probit link can
  be used, depending on the value of \code{outcome_model}.  The model
  assumes that the conditional distribution of \eqn{X}{X} given
  \eqn{(Z,G)}{(Z,G)} is normal with a mean vector that depends on
  \eqn{(Z,G)}{(Z,G)} and a covariance matrix that is assumed not to
  depend on \eqn{(Z,G)}{(Z,G)}.  Both the regression parameters and the
  residual covariance matrix of this conditional distribution are
  estimated.

  All parameters of the model involving \eqn{(Y,X,Z,G)}{(Y,X,Z,G)} are
  estimated using the observed data \eqn{(Y,W,Z,G)}{(Y,W,Z,G)} using
  assumptions and information about the distribution of the measurement
  errors \eqn{U}{U}.  The structure assumed here is that measurement
  errors are independent across units and across dimensions of
  \eqn{X}{X}, and that the conditional distribution of \eqn{U}{U} given
  \eqn{(Y,X,Z,G)} is a normal distribution with mean zero and variance
  \eqn{g(X)}{g(X)}.  The function \eqn{g}{g} must be specified and can
  be constant.  Additional discussions of this class of error functions
  are provided in the references, and details about how information
  about \eqn{g}{g} is conveyed to this function are provided below.
  
  
  \emph{Syntax Details}
  
  Note that this function requires the R2jags package, which in turn
  requires JAGS to be installed on your system.

  The function will check that the only column of \code{Z} that is in
  the span of the columns of the design matrix implied by the grouping
  variable \code{G} is the first column, corresponding to an intercept.
  The effects of \code{G} are parameterized with a sum-to-zero
  constraint, so that the effect of each group is expressed relative to
  the average of all group effects.

  The \code{varfuncs} argument requires the most thinking.  This
  argument is a list with as many elements as there are error-prone
  covariates.  Each element of the list is itself a list providing the
  measurement error information about one of the error-prone covariates,
  call it \code{i} corresponding to \code{W[,i]} to be concrete.  Then
  \code{varfuncs[[i]]} must a list following a particular structure.
  First, \code{varfuncs[[i]]$type} must be a character string taking one
  of three possible values: \code{constant}, \code{piecewise_linear} or
  \code{log_polynomial}.  The case \code{constant} corresponds to the
  case of homoskedastic measurement error where \eqn{g(X)}{g(X)} is
  constant, and the variance of this measurement error must be provided
  in \code{varfuncs[[i]]$vtab}.  The other two cases correspond to the
  case where the conditional measurement error variance \eqn{g(X)}{g(X)}
  is a nontrivial function of \eqn{X}{X}.  In both of these cases,
  \code{varfuncs[[i]]$vtab} must be a matrix or data frame with exactly
  two columns and \code{K} rows, where the first column provides values
  \code{x[1],...,x[K]} of \eqn{X}{X} and the second column provides
  values \code{g(x[1]),...,g(x[K])}.  That is, the function \eqn{g(X)}
  is conveyed via a lookup table.  The value of \code{K} is selected by
  the user. Larger values of \code{K} will make the approximation to
  \eqn{g(X)} more accurate but will cause the model estimation to
  proceed more slowly.  How the values in the lookup table are used to
  approximate \eqn{g(X)} more generally depends whether
  \code{varfuncs[[i]]$type} is \code{piecewise_linear} or
  \code{log_polynomial}.  In the case of \code{piecewise_linear}, the
  values in the lookup table are linearly interpolated.  In the case of
  \code{log_polynomial}, a polynomial of degree
  \code{varfuncs[[i]]$degree} is fitted to the logs of the values of
  \code{g(x[1]),...,g(x[K])}, and the fitted model is used to build a
  smooth approximation to the function \eqn{g(X)}{g(X)}.  The default
  value of \code{varfuncs[[i]]$degree} if it is not specified is 6.  For
  either the piecewise linear or log polynomial approximations, the
  function \code{g(X)}{g(X)} is set to \code{g(x[1])} for values of
  \code{x} smaller than \code{x[1]}, and is set of \code{g(x[K])} for
  values of \code{x} larger than \code{x[K]}.  Diagnostic plots of the
  approximate variance functions saved in PDF file whose location is
  returned by \code{lr_ancova}.  The Examples section provides examples
  that will be helpful in specifying \code{varfuncs}.

  When there are two or more error-prone covariates, the model estimates
  a residual variance/covariance matrix of \eqn{X}{X} given
  \eqn{(Z,G)}{(Z,G)}.  Because the model is fit in a Bayesian framework,
  a prior distribution is required for this matrix.  We are using JAGS
  and specify a prior distribution for the inverse of the residual
  variance/covariance matrix using a Wishart distribution.  The degrees
  of freedom parameter of this distribution is set to one plus
  \code{ncol(W)} to be minimally informative.  The scale matrix of this
  distribution can be set by passing an appropriate matrix via the
  \code{scalemat} argument.  If \code{scalemat} is NULL, the function
  specifies a diagonal scale matrix that attempts to make the prior
  medians of the unknown residual variances approximately equal to the
  residual variances obtained by regressing components of \eqn{W}{W} on
  \eqn{(Z,G)}{(Z,G)}.  See \code{\link{get_bugs_wishart_scalemat}}.
  Such variances will be somewhat inflated due to measurement error in
  \eqn{W}{W} but the prior variance of the Wishart distribution is
  sufficiently large that this lack of alignment should be minimally
  consequential in most applications.  The value of \code{scalemat} used
  in the estimation is returned by the function, and users can start
  with the default and then pass alternative values via the
  \code{scalemat} argument for sensitivity analyses if desired.  }

\value{A object of class \code{rjags}, with additional information
  specific to this context.  The additional information is stored as a
  list called \code{lr_ancova_extras} with the following components:

  \item{model.location}{Path to file containing JAGS model code.}

  \item{plot.location}{Path to file containing diagnostic plots
    regarding the variance functions.}

  \item{group.map}{A dataframe mapping the original group labels in
    \code{G} to integer labels ranging from \code{1} to the number of
    unique elements of \code{G}.  These are useful for mapping the group
    effects reported by JAGS back to the group labels.}

  \item{scalemat}{The value of \code{scalemat} used in the estimation.}

  The parameters used in the JAGS model, and thus named in the model
  object, use naming conventions described here.  The parameters in the
  linear function of \eqn{(X,Z,G)}{(X,Z,G)} that is related to
  \eqn{Y}{Y} are partitioned into \code{betaYXZ} and \code{betaYG}.  In
  applications involving analysis of causal effects of groupings, the
  parameters \code{betaYG} will generally be of most interest.  When
  \code{outcome_model} is \code{normal}, the residual standard deviation
  of \eqn{Y}{Y} given \eqn{(X,Z,G)}{(X,Z,G)} is also estimated and is
  called \code{sdYgivenXZG}.  The regression parameters for the
  conditional distribution of \eqn{X}{X} given \eqn{(Z,G)}{(Z,G)} are
  partitioned as \code{betaXZ} and \code{betaXG}.  The residual
  variance/covariance matrix for \eqn{X}{X} given \eqn{(Z,G)}{(Z,G)} is
  named \code{varXgivenXG}. Additional details on these parameters can
  be found by looking at the JAGS model file whose location is returned
  as noted above.
}

\seealso{
  \code{\link{jags}}, \code{\link{get_bugs_wishart_scalemat}}
}

\references{
  Battauz, M. and Bellio, R. (2011). \dQuote{Structural modeling of
  measurement error in generalized linear models with Rasch measures
  as covariates,} \emph{Psychometrika} 76(1):40-56.
    
  Lockwood J.R. and McCaffrey D.F. (2014). \dQuote{Correcting for test
    score measurement error in ANCOVA models for estimating treatment
    effects,} \emph{Journal of Educational and Behavioral Statistics}
  39(1):22-52.

  Lockwood J.R. and McCaffrey D.F. (2017).
  \dQuote{Simulation-extrapolation with latent heteroskedastic
    variance,} \emph{Psychometrika} 82(3):717-736.

  Plummer, M. (2003). \dQuote{JAGS: A program for analysis of Bayesian
  graphical models using Gibbs sampling.} Proceedings of the 3rd
  International Workshop on Distributed Statistical Computing (DSC
  2003), Vienna, Austria.

  Rabe-Hesketh S., Pickles A. and Skrondal A. (2003). \dQuote{Correcting
  for covariate measurement error in logistic regression using
  nonparametric maximum likelihood estimation,} \emph{Statistical
  Modelling} 3:215-232.
}

\author{
  J.R. Lockwood \email{jrlockwood@ets.org}
}

\examples{
set.seed(3001)
cat("NOTE: this example uses MCMC and takes a little while to run\n")
cat("TODO")
}

\keyword{models}
