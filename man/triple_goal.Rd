\name{triple_goal}
\alias{triple_goal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
triple_goal(s, stop.if.ties = FALSE, quantile.type = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{stop.if.ties}{
%%     ~~Describe \code{stop.if.ties} here~~
}
  \item{quantile.type}{
%%     ~~Describe \code{quantile.type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (s, stop.if.ties = FALSE, quantile.type = 7) 
{
    if (!is.matrix(s)) {
        stop("s must be a matrix")
    }
    if (!is.numeric(s)) {
        stop("s must be numeric")
    }
    if (any(is.na(s))) {
        stop("s must not contain missing values")
    }
    if (!(quantile.type \%in\% 1:9)) {
        stop("quantile.type must be integer between 1 and 9")
    }
    K <- ncol(s)
    theta_pm <- apply(s, 2, mean)
    etadot <- mean(theta_pm)
    theta_psd <- apply(s, 2, sd)
    lambda_k <- theta_psd^2
    var_pm <- var(theta_pm)
    theta_cb <- etadot + (sqrt(1 + (mean(lambda_k)/var_pm)) * 
        (theta_pm - etadot))
    check1 <- mean(theta_cb) - etadot
    check2 <- var(theta_cb) - (mean(lambda_k) + var_pm)
    stopifnot(abs(check1) + abs(check2) < 1e-10)
    rbar <- apply(t(apply(s, 1, rank, ties.method = "average")), 
        2, mean)
    if (any(diff(sort(rbar)) < .Machine$double.eps) && stop.if.ties) {
        stop("posterior mean ranks not unique")
    }
    rhat <- rank(rbar, ties.method = "random")
    stopifnot(all(sort(rhat) == 1:K))
    theta_gr <- quantile(c(s), probs = (2 * rhat - 1)/(2 * K), 
        names = FALSE, type = quantile.type)
    e <- data.frame(index = 1:K, theta_pm = theta_pm, theta_psd = theta_psd, 
        theta_cb = theta_cb, theta_gr = theta_gr, rbar = rbar, 
        rhat = rhat)
    if (!is.null(colnames(s))) {
        rownames(e) <- colnames(s)
    }
    else {
        rownames(e) <- 1:nrow(e)
    }
    return(e)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
