#' Density Estimation based on Distribution Element Trees
#'
#' The function \code{det.query} evaluates probability densities at the query points \code{x} based on a distribution element tree (DET). The latter is calculable with \code{\link{det.construct}} based on available data.
#'
#' @param det distribution element tree object resulting from \code{\link{det.construct}} based on data with \code{d} components or dimensions.
#' @param x matrix containing \code{n} query points (columns) with \code{d} components or dimensions (rows).
#' @param cores for large query-point sets, \code{cores > 1} allows for parallel tree query. The default is \code{cores = 0}, which allocates half of the available cores (see \code{\link[parallel]{detectCores}}). \code{cores = 1} corresponds to serial tree query.
#'
#' @return A vector containing the probability density at the query points \code{x} is returned.
#' @export
#'
#' @examples
#' ## 1d example
#' require(stats); require(graphics)
#' # DET generation based on Gaussian data
#' det <- det.construct(matrix(rnorm(1e5,2,3), nrow = 1), cores = 2)
#' # density evaluation based on DET at equidistant query points
#' x <- matrix(seq(-10,14,0.01), nrow = 1); p <- det.query(det, x, cores = 1)
#' # compare DET estimate (black) against Gaussian reference (red)
#' plot(x, p, type = "l", col = "black")
#' lines(x, dnorm(x,2,3), col = "red")
#'
#' ## 2d example
#' require(stats); require(graphics)
#' # mean and covariance of Gaussian, data generation
#' mu <- c(3,5); C <- matrix(c(4.0,-2.28,-2.28,1.44), nrow = 2)
#' A <- eigen(C); B <- diag(A$values); A <- A$vectors
#' x <- matrix(rnorm(2e4), nrow = 2)
#' x <- t(A %*% (sqrt(B) %*% x) + mu %*% t(rep(1,dim(x)[2])))
#' # bounds and resolution of x1-x2 query grid
#' lb <- c(-5,0); ub <- c(11,10); np <- c(320,200)
#' x1 <- lb[1] + (ub[1]-lb[1])*((1:np[1])-0.5)/np[1]
#' x2 <- lb[2] + (ub[2]-lb[2])*((1:np[2])-0.5)/np[2]
#' xp <- rbind(rep(x1,np[2]), rep(x2,each = np[1])) # grid points
#' # plotting
#' split.screen(c(2, 2)); screen(1)
#' plot(x, type = "p", pch = ".", asp = 1, main = "data")
#' # DET estimator
#' det <- det.construct(t(x), cores = 1)
#' yd <- matrix(det.query(det, xp, cores = 2), nrow = np[1])
#' screen(2)
#' image(list(x = x1, y = x2, z = yd), asp = 1,
#'       col = grDevices::gray((100:0)/100), main = "det")
#' # Gaussian density for comparison
#' yr <- yr <- exp(-1/2 * colSums(
#'    (t(solve(C)) %*% (xp - mu%*%t(rep(1,dim(xp)[2])))) *
#'                     (xp - mu%*%t(rep(1,dim(xp)[2]))))
#'                               ) / sqrt((2*pi)^2*det(C))
#' yr <- matrix(yr, nrow = np[1])
#' screen(3)
#' image(list(x = x1, y = x2, z = yr), asp = 1,
#'       col = grDevices::gray((100:0)/100), main = "reference")
det.query <- function(det, x, cores = 0) {
   # check data
   errorstrg <- sprintf(
      "invalid x, expecting d x n matrix with d = %i dimensions and n query points",
      length(det$lb))
   if (is.matrix(x)) {
      if (dim(x)[2] == length(det$lb)) {stop(errorstrg)}
   } else {stop(errorstrg)}
   n <- dim(x)[2] # number of points
   # transform and normalize query points
   x <- (det$A%*%(x-det$mu%*%t(rep(1,n))) - det$lb%*%t(rep(1,n))) /
      ((det$ub - det$lb)%*%t(rep(1,n))) # x in [0,1]
   # setup cluster for parallel processing
   if (cores <= 0) {cores <- max(1,ceiling(parallel::detectCores()/2))}
   clster <- parallel::makeCluster(cores)
   # loop over points in x
   x <- split(x, col(x)) # matrix -> list of column vectors
   p <- parallel::parLapply(clster, x, function(xk) {
      # check if point xk is outside of bounds
      if (any(xk < 0) | any(xk > 1)) {pk <- 0} # outside det bounds -> p = 0
      else { # inside det bounds -> find de
         ind <- 1 # start from root de
         while (!is.na(det$tree[ind,2])) { # loop until leaf or final de is found
            dimens <- det$sd[ind]; pos <- det$sp[ind] # split dimension & position
            # xk in first or second child de?
            if (xk[dimens] <= pos) { # first child
               ind <- det$tree[ind,2] # de index
               xk[dimens] <- xk[dimens] / pos # renormalize component
            } else { # second child
               ind <- det$tree[ind,3] # de index
               xk[dimens] <- (xk[dimens]-pos) / (1-pos) # renormalize component
            }
         }
         # determine probability density
         pk <- det$p[ind]
         if (pk > 0) {pk <- pk * prod(((xk-1/2)*det$theta[[ind]]+1))}
      }
      return(pk)
   })
   parallel::stopCluster(clster)
   # rescale
   p <- unlist(p) / prod(det$ub - det$lb)
   return(p)
}
