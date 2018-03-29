#' Density Estimation Based on Distribution Element Trees
#'
#' The function \code{det.query} evaluates probability densities at the query points \code{x} based on a distribution element tree (DET). The latter is calculable with \code{\link{det.construct}} based on available data.
#'
#' @param det distribution element tree object resulting from \code{\link{det.construct}} based on data with \code{d} components or dimensions.
#' @param x matrix containing \code{n} query points (columns) with \code{d} components or dimensions (rows).
#' @param cores for large query-point sets, \code{cores > 1} allows for parallel tree query using the indicated number of cores. \code{cores = Inf} allocates half of the available cores (see \code{\link[parallel]{detectCores}}). The default is \code{cores = 1} corresponding to serial tree query.
#'
#' @return A vector containing the probability density at the query points \code{x} is returned.
#' @export
#'
#' @examples
#' ## 1d example
#' require(stats); require(graphics)
#' # DET generation based on Gaussian/uniform data
#' det <- det.construct(t(c(rnorm(1e5,2,3),runif(1e5)-3)))
#' # density evaluation based on DET at equidistant query points
#' x <- t(seq(-10,14,0.01)); p <- det.query(det, x)
#' # compare DET estimate (black) against Gaussian/uniform reference (red)
#' plot(x, p, type = "l", col = "black")
#' lines(x, (dnorm(x,2,3)+dunif(x+3))/2, col = "red")
#'
#' ## 2d example
#' require(stats); require(graphics)
#' # mean and covariance of Gaussian, data generation
#' mu <- c(3,5); C <- matrix(c(4.0,-2.28,-2.28,1.44), nrow = 2)
#' A <- eigen(C); B <- diag(A$values); A <- A$vectors
#' x <- matrix(rnorm(2e4), nrow = 2)
#' x <- t(A %*% (sqrt(B) %*% x) + mu %*% t(rep(1,ncol(x))))
#' # bounds and resolution of x1-x2 query grid
#' lb <- c(-5,0); ub <- c(11,10); np <- c(320,200)
#' x1 <- lb[1] + (ub[1]-lb[1])*((1:np[1])-0.5)/np[1]
#' x2 <- lb[2] + (ub[2]-lb[2])*((1:np[2])-0.5)/np[2]
#' xp <- rbind(rep(x1,np[2]), rep(x2,each = np[1])) # grid points
#' # plotting
#' split.screen(c(2, 2)); screen(1)
#' plot(x, type = "p", pch = ".", asp = 1, main = "data")
#' # DET estimator
#' det <- det.construct(t(x))
#' yd <- matrix(det.query(det, xp), nrow = np[1])
#' screen(2)
#' image(list(x = x1, y = x2, z = yd), asp = 1,
#'       col = grDevices::gray((100:0)/100), main = "det")
#' # Gaussian density for comparison
#' yr <- yr <- exp(-1/2 * colSums(
#'    (t(solve(C)) %*% (xp - mu%*%t(rep(1,ncol(xp))))) *
#'                     (xp - mu%*%t(rep(1,ncol(xp)))))
#'                               ) / sqrt((2*pi)^2*det(C))
#' yr <- matrix(yr, nrow = np[1])
#' screen(3)
#' image(list(x = x1, y = x2, z = yr), asp = 1,
#'       col = grDevices::gray((100:0)/100), main = "reference")
det.query <- function(det, x, cores = 1) {
   # check data
   errorstrg <- sprintf(
      "invalid x, expecting d x n matrix with d = %i dimensions and n query points",
      length(det$lb))
   if (is.matrix(x)) {
      if (ncol(x) == length(det$lb)) {stop(errorstrg)}
   } else {stop(errorstrg)}
   n <- ncol(x) # number of points
   # transform and normalize query points
   x <- (det$A%*%(x-det$mu%*%t(rep(1,n))) - det$lb%*%t(rep(1,n))) /
      ((det$ub - det$lb)%*%t(rep(1,n))) # x in [0,1]
   x <- split(x, col(x)) # for list-apply, matrix -> list of column vectors
   # determine density based on det estimator at point x
   ddet <- function(x) {
      # check if point x is outside of bounds
      if (any(x < 0) | any(x > 1)) {p <- 0} # outside det bounds -> p = 0
      else { # inside det bounds -> find de
         ind <- 1 # start from root de
         while (!is.na(det$tree[ind,2])) { # loop until leaf or final de is found
            dimens <- det$sd[ind]; pos <- det$sp[ind] # split dimension & position
            # x in first or second child de?
            if (x[dimens] <= pos) { # first child
               ind <- det$tree[ind,2] # de index
               x[dimens] <- x[dimens] / pos # renormalize component
            } else { # second child
               ind <- det$tree[ind,3] # de index
               x[dimens] <- (x[dimens]-pos) / (1-pos) # renormalize component
            }
         }
         # determine probability density
         p <- det$p[ind]
         if (p > 0) {p <- p * prod(((x-1/2)*det$theta[[ind]]+1))}
      }
      return(p)
   }
   # loop over points in x
   if (cores > 1) { # parallel processing
      if (is.infinite(cores)) {cores <- max(1,ceiling(parallel::detectCores()/2))}
      clster <- parallel::makeCluster(cores)
      p <- parallel::parLapply(clster, x, ddet)
      parallel::stopCluster(clster)
   } else { # serial
      p <- lapply(x, ddet)
   }
   # rescale
   p <- unlist(p) / prod(det$ub - det$lb)
   return(p)
}
