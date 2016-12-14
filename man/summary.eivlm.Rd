\name{summary.eivlm}
\alias{summary.eivlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.eivlm(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{correlation}{
%%     ~~Describe \code{correlation} here~~
}
  \item{symbolic.cor}{
%%     ~~Describe \code{symbolic.cor} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (object, correlation = FALSE, symbolic.cor = FALSE, 
    ...) 
{
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
        class(ans) <- "summary.eivlm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0L, 4L)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
            "Std. Error", "t value", "Pr(>|t|)"))
        ans$latent_resvar <- rss/rdf
        ans$r.squared <- 0
        return(ans)
    }
    if (is.null(z$terms)) 
        stop("invalid 'eivlm' object:  no 'terms' component")
    if (!inherits(object, "eivlm")) 
        warning("calling summary.eivlm(<fake-eivlm-object>) ...")
    if (is.null(Qr <- z$unadj_qr)) {
        stop("eivlm object does not have a proper 'qr' component for unadjusted model.\n Rank zero or should not have used eivreg(.., qr=FALSE).")
    }
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual) 
        warning("residual degrees of freedom in object suggest this is not an \"eivlm\" fit")
    r <- z$unadj_residuals
    f <- z$unadj_fitted.values
    w <- z$weights
    if (is.null(w)) {
        unadj_mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        unadj_rss <- sum(r^2)
    }
    else {
        unadj_mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        unadj_rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    unadj_resvar <- unadj_rss/rdf
    if (is.finite(unadj_resvar) && unadj_resvar < (mean(f)^2 + 
        var(f)) * 1e-30) 
        warning("essentially perfect fit for unadjusted model: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * unadj_resvar)
    est <- z$unadj_coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$unadj_residuals <- r
    ans$unadj_coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
        rdf, lower.tail = FALSE))
    dimnames(ans$unadj_coefficients) <- list(names(z$unadj_coefficients)[Qr$pivot[p1]], 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$unadj_sigma <- sqrt(unadj_resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept")) 
            1L
        else 0L
        ans$unadj_r.squared <- unadj_mss/(unadj_mss + unadj_rss)
        ans$unadj_adj.r.squared <- 1 - (1 - ans$unadj_r.squared) * 
            ((n - df.int)/rdf)
        ans$unadj_fstatistic <- c(value = (unadj_mss/(p - df.int))/unadj_resvar, 
            numdf = p - df.int, dendf = rdf)
    }
    else {
        ans$unadj_r.squared <- ans$unadj_adj.r.squared <- 0
    }
    ans$unadj_cov.unscaled <- R
    dimnames(ans$unadj_cov.unscaled) <- dimnames(ans$unadj_coefficients)[c(1, 
        1)]
    if (correlation) {
        ans$unadj_correlation <- (R * unadj_resvar)/outer(se, 
            se)
        dimnames(ans$unadj_correlation) <- dimnames(ans$unadj_cov.unscaled)
        ans$unadj_symbolic.cor <- symbolic.cor
    }
    est <- z$coefficients
    passthrough <- c("residuals", "fitted.values", "N", "latent_resvar", 
        "vcov", "relnames")
    ans[passthrough] <- z[passthrough]
    if (!is.null(z$reliability)) {
        ans$reliability <- z$reliability
    }
    else {
        ans$Sigma_error <- z$Sigma_error
    }
    ans$latent_R2_dfadj <- (z$varYXZ[1, 1] - z$latent_resvar)/z$varYXZ[1, 
        1]
    ans$latent_R2 <- (z$varYXZ[1, 1] - ((z$latent_resvar) * rdf/(z$N)))/z$varYXZ[1, 
        1]
    if (!is.null(z$cluster_varname)) {
        passthrough <- c("cluster_varname", "cluster_values", 
            "cluster_num")
        ans[passthrough] <- z[passthrough]
    }
    if (!is.null(ans$vcov)) {
        se <- sqrt(diag(ans$vcov))
        tval <- est/se
        ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
            rdf, lower.tail = FALSE))
        dimnames(ans$coefficients) <- list(rownames(ans$vcov), 
            c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        if (correlation) {
            ans$correlation <- ans$vcov/outer(se, se)
            dimnames(ans$correlation) <- dimnames(ans$vcov)
        }
    }
    else {
        ans$coefficients <- cbind(est, NA, NA, NA)
        dimnames(ans$coefficients) <- list(names(est), c("Estimate", 
            "Std. Error", "t value", "Pr(>|t|)"))
        if (correlation) {
            ans$correlation <- matrix(0, nrow = length(est), 
                ncol = length(est))
            is.na(ans$correlation) <- TRUE
            dimnames(ans$correlation) <- list(names(est), names(est))
        }
    }
    ans$symbolic.cor <- symbolic.cor
    if (!is.null(z$na.action)) 
        ans$na.action <- z$na.action
    class(ans) <- "summary.eivlm"
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
