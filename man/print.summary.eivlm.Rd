\name{print.summary.eivlm}
\alias{print.summary.eivlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print method for \code{summary.eivlm} objects}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary.eivlm(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{symbolic.cor}{
%%     ~~Describe \code{symbolic.cor} here~~
}
  \item{signif.stars}{
%%     ~~Describe \code{signif.stars} here~~
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
function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (!is.null(x$reliability)) {
        cat("Reliability:\n")
        print(x$reliability, digits = digits)
    }
    else {
        cat("Error Covariance Matrix\n")
        print(x$Sigma_error, digits = digits)
    }
    cat("\n")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if (!is.null(x$weights) && diff(range(x$weights))) 
        "Weighted ", "Residuals:\n", sep = "")
    if (rdf > 5L) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2L) {
            structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                dimnames(resid)[[2L]]))
        }
        else {
            zz <- zapsmall(quantile(resid), digits + 1L)
            structure(zz, names = nam)
        }
        print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
        print(resid, digits = digits, ...)
    }
    else {
        cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
        cat("\n")
    }
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    }
    else {
        if (nsingular <- df[3L] - df[1L]) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("Number of observations used:", x$N, "\n")
    if (!is.null(x$cluster_varname)) {
        cat(paste("Std. Errors adjusted for ", x$cluster_num, 
            " clusters in ", x$cluster_varname, "\n", sep = ""))
    }
    cat("Latent residual standard deviation:", format(signif(sqrt(x$latent_resvar), 
        digits)), "\n")
    cat("Latent R-squared: ", format(signif(x$latent_R2, digits)), 
        ", (df-adjusted: ", format(signif(x$latent_R2_dfadj, 
            digits)), ")\n", sep = "")
    cat("\nEIV-Adjusted vs Unadjusted Coefficients:\n")
    tmp <- cbind(x$coefficients[, 1, drop = FALSE], x$unadj_coefficients[, 
        1])
    colnames(tmp) <- c("Adjusted", "Unadjusted")
    print(tmp, digits = digits)
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
