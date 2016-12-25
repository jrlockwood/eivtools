deconv_npmle <- function(W, csem, xmin = -5, xmax = +5, ngrid = 2000, lambda = 0.0005, lltol = 1e-7, psmall = 0.00005, discrete = FALSE, quietly = FALSE){
    
    stopifnot( all(!is.na(W)) && is.function(csem) && (xmin < xmax) && (lambda > 0) && (lltol > 0) )
    if(discrete){
        stop("discrete option not currently implemented")
    }
    if(min(W) < xmin){
        stop("xmin exceeds smallest data value, consider decreasing xmin")
    }
    if(max(W) > xmax){
        stop("largest data value exceeds xmax; consider increasing xmax")
    }
    
    ## NOTE: we collapse over repeat W values when possible, using "counts" to contribute to likelihood
    .W <- sort(unique(W))
    nW <- length(.W)
    .counts <- sapply(.W, function(x){ sum(W == x) })
    stopifnot(sum(.counts) == length(W))
    W <- .W
    
    ## functions to map probabilities to reduced-dimension unconstrained scale, and inverse
    theta_to_p <- function(theta){
        if(length(theta) == 0){
            return(1)
        } else {
            p <- exp(theta) / (1 + sum(exp(theta)))
            return(c(p, 1-sum(p)))
        }
    }
    
    p_to_theta <- function(p){
        stopifnot( all(p > 0) && (abs((1 - sum(p))) < 1e-9) )
        if(length(p) == 1){
            return(numeric(0))
        } else {
            last <- length(p)
            return(log(p[1:(last-1)] * (1 + ((1-p[last])/p[last])) ))
        }
    }
  
    ## CHECKS:
    ## p <- c(0.1, 0.2, 0.4, 0.3); print(ma(p - theta_to_p(p_to_theta(p))))
    ## theta <- c(-3.2, 1.0, -1.0); print(ma(theta - p_to_theta(theta_to_p(theta))))
    ## print(ma(1 - theta_to_p(p_to_theta(1))))
    ## p_to_theta(theta_to_p(numeric(0)))
    
    ## build (nW x ngrid) matrix of conditional densities
    grid      <- seq(from=xmin, to=xmax, length=ngrid)
    grid.csem <- csem(grid)
    ## if(!discrete){
    fwx <- dnorm( matrix(W,ncol=ngrid,nrow=nW), mean=matrix(grid,ncol=ngrid,nrow=nW,byrow=T), sd=matrix(grid.csem,ncol=ngrid,nrow=nW,byrow=T))
    ##  } else { ## uses pxgu - fwx[i,j] = p(W=W[i] | X = grid[j])
    ## tmp <- lapply(grid, function(u){ pxgu(u, csem(u), .W)})
    ## stopifnot(all(sapply(tmp, function(x){ all(x[,1] == .W) })))
    ## fwx <- do.call("cbind", lapply(tmp, function(x){ x[,2] }))
    ## }
  
    ## negative log likelihood for a set of probabilities given ".inds" which are
    ## indices of "grid" and which are continually updated
    ## NOTE: updated to use counts
    negll <- function(theta){
        -sum(.counts * log(as.vector(fwx[,.inds] %*% theta_to_p(theta))))
    }
    
    ## find initial best grid point.  NOTE this will often want the X with the
    ## biggest CSEM because a single point is inconsistent with most W unless the
    ## CSEM is large. the exception is if we pick ugrid to be very large, then it
    ## will find interior point.  however even if it picks an extreme starting
    ## point, that point will be dropped later if it is too far in the tails.
    ll <- apply(fwx, 2, function(x){ sum(.counts * log(x)) })
    .inds <- which.max(ll)
    .probs   <- 1
    .eligible <- setdiff(1:ngrid, .inds)
    ll.current <- ll[.inds]
    
    .history <- list()
    .history[[1]] <- data.frame(x = grid[.inds], csem = grid.csem[.inds], p = 1, ll = ll.current, ex = grid[.inds], varx = 0)

    ## now add grid points per the algorithm until there is no improvement
    done <- FALSE  
    while(!done){
        ## evaluate each eligible point with a weight "lambda"
        part.static   <- as.vector(fwx[,.inds,drop=FALSE] %*% ((1 - lambda)*.probs))
        part.total   <- matrix(part.static, ncol=length(.eligible), nrow=nW) + (fwx[,.eligible] * lambda)
        ll.candidate <- apply(part.total, 2, function(x){ sum(.counts * log(x)) })
        if(all(ll.candidate - ll.current < lltol)){
            done <- TRUE
        } else {
            w <- which.max(ll.candidate - ll.current)
            .inds <- c(.inds, .eligible[w])
            .eligible <- setdiff(.eligible, .eligible[w])
            
            ## set starting value: mass 0.05 at new value, normalized p on existing values
            o <- optim(p_to_theta(c(0.95*(.probs/sum(.probs)), 0.05)), negll, method  = "BFGS", control = list(maxit=1000, trace=5*(1 - as.numeric(quietly)), REPORT = 1))
            .probs <- theta_to_p(o$par)
            ll.current <- -o$value
            
            ## sometimes we might have picked up an early grid point
            ## but now it has virtually no probability, so dump it
            w <- which(.probs < psmall)
            if(length(w) > 0){
                .probs <- .probs[-w]
                .probs <- .probs/sum(.probs)
                .inds  <- .inds[-w]
            }
            
            ## reorder the grid if need be
            o <- order(.inds)
            .inds <- .inds[o]
            .probs <- .probs[o]
            
            ## summarize the state
            .history[[length(.history)+1]] <-
                data.frame(x = grid[.inds],
                           csem = grid.csem[.inds],
                           p = .probs,
                           ll = ll.current,
                           ex = sum(grid[.inds] * .probs),
                           varx = sum(grid[.inds]^2 * .probs) - (sum(grid[.inds] * .probs))^2)
            if(!quietly){
                print(.history[[length(.history)]])
                cat("\n\n\n")      
            }
        }
    }
    return(list(.history = .history, px = .history[[length(.history)]]))
}
