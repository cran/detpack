#' Bootstrapping from Distribution Element Trees
#'
#' Smooth bootstrapping or generation of (un)conditional random vectors based on an existing distribution element tree (DET).
#'
#' @param n number of samples to generate.
#' @param det distribution element tree object resulting from \code{\link{det.construct}}.
#' @param xc vector with conditioning values of probability-space components listed in \code{dc}. If empty (default), unconditional samples are generated.
#' @param dc integer vector with indices of conditioning components corresponding to \code{xc}. If empty (default), unconditional samples are generated.
#' @param cores for large \code{n}, \code{cores > 1} allows for parallel bootstrapping. The default is \code{cores = 0}, which allocates half of the available cores (see \code{\link[parallel]{detectCores}}). \code{cores = 1} corresponds to serial bootstrapping.
#'
#' @return A matrix containing \code{n} random vectors (columns) with \code{d} components or dimensions (rows) is returned. \code{d} is equal to the dimensionality of the underlying \code{det} object.
#' @export
#'
#' @examples
#' ## 2d example
#' require(stats); require(graphics)
#' # data from uniform distribution on a wedge
#' x <- matrix(runif(2e4), ncol = 2); x <- x[x[,2]<x[,1],]
#' x2c <- 0.75 # conditioning component
#' # data and conditioning line
#' split.screen(c(2, 1)); screen(1)
#' plot(x, type = "p", pch = ".", asp = 1)
#' lines(c(0,1), x2c*c(1,1), col = "red")
#' # DET construction and bootstrapping
#' det <- det.construct(t(x), mode = 1, lb = 0, ub = 0, cores = 1) # const. de's, no pre-white
#' y <- det.rnd(1e3, det, xc = x2c, dc = 2, cores = 1) # conditional bootstrap'g
#' # compare generated data (black) with exact cond. distribution (red)
#' screen(2); det1(y[1,], mode = 1, cores = 2)
#' lines(c(0,x2c,x2c,1,1),c(0,0,1/(1-x2c),1/(1-x2c),0), col = "red")
#'
#' grDevices::dev.off() # erase plot
#'
#' ## example 2d unconditional
#' require(stats); require(graphics)
#' x <- matrix(runif(1e4), ncol = 2); x <- x[x[,2]<x[,1],] # uniform wedge
#' det <- det.construct(t(x), mode = 1, lb = 0, ub = 0, cores = 1) # no pre-white
#' y <- t(det.rnd(nrow(x), det, cores = 1)) # smooth bootstrapping
#' split.screen(c(2, 1))
#' screen(1); plot(x, type = "p", pch = ".", asp = 1, main = "original")
#' screen(2); plot(y, type = "p", pch = ".", asp = 1, main = "bootstrapped")
#'
#' grDevices::dev.off() # erase plot
#'
#' ## example 3d
#' require(stats); require(graphics)
#' # mean and covariance of Gaussian, data generation
#' mu <- c(1,3,2); C <- matrix(c(25,7.5,1.75,7.5,7,1.35,1.75,1.35,0.43), nrow = 3)
#' A <- eigen(C); B <- diag(A$values); A <- A$vectors
#' x <- matrix(rnorm(24e3), nrow = 3)
#' x <- A %*% (sqrt(B) %*% x) + mu %*% t(rep(1,ncol(x)))
#' lbl <- "x1 | x2 = 7 & x3 = 2.5"
#' pairs(t(x), labels = c("x1","x2","x3"), pch = ".", main = lbl)
#' # bootstrapping conditional on x2 and x3
#' det <- det.construct(x, lb = 0, ub = 0, cores = 2)
#' xc <- c(2.5,7); d <- c(3,2) # conditional on x2 = 7 & x3 = 2.5
#' y <- det.rnd(1e3, det, xc, d, cores = 2)
#' det1(y[1,], mode = 1, main = lbl, cores = 2)
#' # compare with exact conditional density
#' Cm1 <- solve(C); var1 <- det(C)/det(C[2:3,2:3]) # conditional variance
#' mu1 <- mu[1] + var1*((mu[2]-xc[d==2])*Cm1[1,2]+(mu[3]-xc[d==3])*Cm1[1,3]) # cond. mean
#' x1 <- mu1 + seq(-50,50)/50 * 5*sqrt(var1) # x1-axis grid points
#' lines(x1, dnorm(x1,mu1,sqrt(var1)), col = "red")
det.rnd <- function(n, det, xc = vector("numeric", length = 0), dc = vector("numeric", length = 0), cores = 0) {
   ds <- length(det$lb) # number of sample-space dimensions
   # collect det leafs cut by condit'g planes p
   if (length(xc) == 0) { # unconditional case
      leafs <- seq(nrow(det$tree))[is.na(det$tree[,2])] # list with indices of leaf elements
   } else { # conditional case
      # check input
      if (!(all(det$A == diag(ds)) & all(det$mu == 0))) {
         stop("pre-whitening not supported (det$A != I, use det.construct with lb = ub = 0)")
      }
      dc <- unique(dc) # remove possible duplicate dimensions
      if (!((length(xc) == length(dc)) & all(ds >= dc) & all(dc >= 1))) {
         stop("invalid conditioning inputs xc or dc")
      }
      # extract cut tree leafs
      leafs <- det.cut(det, xc, dc) # list with indices of cut elements
   }
   if (length(leafs) == 0) {return(NULL)} # condit'g planes outside sample space
   # extract leaf properties
   m <- length(leafs) # number of cut leafs
   p <- rep(NA,m) # de probability density
   theta <- rep(list(NA),m) # de parameters
   lb <- matrix(rep(NA,ds*m), nrow = ds) # de lower bound
   size <- lb # de size
   for (k in 1:m) { # loop over cut leafs
      dek <- det.de(det, leafs[k])
      p[k] <- dek$p; theta[[k]] <- dek$theta
      lb[,k] <- dek$lb; size[,k] <- dek$size
   }
   # determine (conditional) leaf probabilites
   sl <- matrix(rep(0,ds*m), nrow = ds) # slope parameter
   for (k in 1:m) {if (!is.na(theta[k])) {sl[,k] <- theta[[k]]}}
   p <- p * apply(size, 2, prod) *
      apply((((xc%*%t(rep(1,m))-lb[dc,])/size[dc,] - 1/2)*sl[dc,] + 1) / size[dc,], 2, prod)
   # setup cluster for parallel processing
   if (cores <= 0) {cores <- max(1,ceiling(parallel::detectCores()/2))}
   clster <- parallel::makeCluster(cores)
   # generate random samples
   p <- cumsum(p)
   if (p[m] == 0) {return(NULL)} # zero probability
   p <- c(0,p/p[m]) # normalize probabilities
   rl <- stats::runif(n*(ds-length(dc)+1))
   rl <- matrix(rl, ncol = n); rl[1,rl[1,] == 0] <- 1
   x <- parallel::parCapply(clster, rl, function(r) {
      x <- rep(0,ds)
      ind <- seq(m)[(p[1:m] < r[1]) & (r[1] <= p[2:(m+1)])]
      # randomly draw sample inside selected leaf with inverse marginal cdfs
      k <- 1 # counter of unconditional, random directions
      for (dimens in 1:ds) { # loop over sample-space dimensions
         j <- seq(along = dc)[dc == dimens] # index of dimens in condit'g list
         if (length(j) == 0) { # dimens is not condit'g direction
            x0 <- lb[dimens,ind]; x1 <- lb[dimens,ind]+size[dimens,ind]
            if (is.na(theta[ind])) {sl <- 0} else {sl <- theta[[ind]][dimens]}
            k <- k+1 # increment counter
            if (sl == 0) { # constant de with uniform pdf
               x[dimens] <- x0 + (x1-x0)*r[k]
            } else { # linear de with linear pdf
               # generate sample from inverse cdf
               a <- sqrt((4 + sl*(sl+8*r[k]-4)) / (x1-x0)^2)
               x[dimens] <- 1/(2*sl) * (x0^2*a + x0*(sl+2-2*a*x1) + x1*(sl-2+a*x1))
            }
         } else { # dimens is condit'g direction
            x[dimens] <- xc[j] # assign conditional value
         }
      }
      return(x)
   })
   parallel::stopCluster(clster)
   x <- matrix(x, nrow = ds)
   # in unconditional case, transform samples back
   if (length(xc) == 0) {x <- t(det$A)%*%x + det$mu%*%t(rep(1,n))}
   return(x)
}

#' Identify Tree Leafs Intersected by Condition(s)
#'
#' Identify distribution element tree (DET) leafs that are cut by conditions. The latter are defined in terms of positions \code{xc} along probability-space components with indices \code{dc}.
#'
#' @param det distribution element tree object resulting from \code{\link{det.construct}}.
#' @param xc vector with conditioning values of probability-space components listed in \code{dc}.
#' @param dc integer vector with indices of conditioning components corresponding to \code{xc}.
#'
#' @return A vector containing the leaf indices that are cut by conditions \code{xc} of components \code{dc} is returned. If no leafs are found, the return vector has length 0.
#' @export
#'
#' @examples
#' # DET based on Gaussian data
#' require(stats); require(graphics)
#' n <- 8e4; x <- rnorm(n)
#' x <- matrix(c(x, x+rnorm(n,0,0.2)), ncol = 2)
#' det <- det.construct(t(x), lb = 0, ub = 0, cores = 2) # no pre-whitening
#' plot(x, type = "p", pch = ".", asp = 1)
#' # leaf elements that are cut by x1 = 2
#' leafs <- det.cut(det, xc = 2, dc = 1) # condition x1 = 2
#' # draw probability space (black) with cut leaf elements (red)
#' rect(det$lb[1], det$lb[2], det$ub[1], det$ub[2], border = "black")
#' for (k in 1:length(leafs)) {
#'    p <- det.de(det, leafs[k])$lb; w <- det.de(det, leafs[k])$size
#'    rect(p[1],p[2],p[1]+w[1],p[2]+w[2], border = "red")
#' }
#' # leafs cut by two conditions x1 = -3, x2 = -2 (blue)
#' leafs <- det.cut(det, xc = c(-2,-3), dc = c(2,1))
#' p <- det.de(det, leafs[1])$lb; w <- det.de(det, leafs[1])$size
#' rect(p[1],p[2],p[1]+w[1],p[2]+w[2], border = "blue")
det.cut <- function(det, xc, dc) {
   m <- 0; leafs <- vector("numeric", length = m) # counter and list of cut det leafs
   xc <- (xc - det$lb[dc]) / (det$ub[dc] - det$lb[dc]) # xc in [0,1]
   # no leafs are cut with condit'g planes outside of sample space
   if (any((xc < 0) | (1 < xc))) {return(leafs)}
   # initialize list of interim elements for search of cut final elements
   ind <- rep(NA,nrow(det$tree)); n <- 0 # list and counter of interim elements
   pr <- matrix(rep(NA,length(xc)*nrow(det$tree)), nrow = length(xc)) # relative condit'g plane positions
   ind[1] <- 1; n <- 1; pr[,1] <- xc # start from root de
   # search det level-by-level for leafs that are cut by condit'g planes
   leafs <- rep(NA,nrow(det$tree))
   while (n > 0) { # loop until all interim elements were searched
      nc <- n # current number of interim elements to be checked
      k <- 1 # counter for ind list
      while (k <= nc) { # loop over interim elements
         # is element a tree leaf or final element?
         if (is.na(det$tree[ind[k],2])) { # yes, final
            # add to list of cut det leafs
            m <- m+1; leafs[m] <- ind[k]
            # remove leaf from list of leafs to be checked
            if (n > k) {ind[k:(n-1)] <- ind[(k+1):n]; pr[,k:(n-1)] <- pr[,(k+1):n]}
            nc <- nc-1; n <- n-1
         } else { # no, just interim
            # check children of interim element
            dimens <- det$sd[ind[k]]; pos <- det$sp[ind[k]] # split dimension & position
            if (any(dc == dimens)) { # split along one of the condit'g planes?
               # replace interim element by child cut by condit'g plane
               if (pr[dc == dimens,k] <= pos) { # first child
                  ind[k] <- det$tree[ind[k],2]
                  pr[dc == dimens,k] <- pr[dc == dimens,k]/pos
               } else { # second child
                  ind[k] <- det$tree[ind[k],3]
                  pr[dc == dimens,k] <- (pr[dc == dimens,k]-pos)/(1-pos)
               }
            } else { # split does not interfere with condit'g -> proceed with both children
               n <- n+1; ind[n] <- det$tree[ind[k],3]; pr[,n] <- pr[,k] # add de's second child
               ind[k] <- det$tree[ind[k],2] # replace interim de by de's first child
            }
            k <- k+1
         }
      }
   }
   return(leafs[1:m])
}
