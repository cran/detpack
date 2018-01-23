#' Density Estimation for Univariate Data Based on Distribution Element Trees
#'
#' One-dimensional piecewise linear or constant probability density estimator based on distribution element trees (DETs).
#'
#' @param x vector with data points.
#' @param mode order of distribution elements applied, default is \code{mode = 2}. Use \code{+/-1} for constant or \code{+/-2} for linear elements. \code{mode > 0} and \code{mode < 0} lead to equal-size and -score splits, respectively, in the element-refinement process.
#' @param bounds \code{c(lb,ub)}, where \code{lb} and \code{ub} are lower and upper bounds of the probability space. If both bounds are set to 0 (default), the bounds are determined based on the data \code{x}.
#' @param alpha significance level for goodness-of-fit testing in element refinement or splitting process. Default is \code{alpha = 1.0e-3}.
#' @param main an overall plot title, see \code{\link[graphics]{title}}.
#' @param cores number of cores for parallel tree construction. Default is \code{0} for automatic allocation, see \code{\link{det.construct}}.
#'
#' @export
#'
#' @examples
#' require(stats)
#' det1(rnorm(1e4), mode = 2, cores = 1)
#' det1(runif(1e4), mode = -1, bounds = c(-1,1.5), cores = 2)
#' det1(rbeta(1e4, shape1 = 1.05, shape2 = 0.8), mode = -1, bounds = c(-0.1,1.1), cores = 1)
#' det1(rbeta(1e4, shape1 = 1.05, shape2 = 0.8), mode = -2, bounds = c(-0.1,1.1), cores = 1)
det1 <- function(x, mode = 2, bounds = c(0,0), alpha = 1.0e-3, main = NULL, cores = 0) {
   # check data
   if (!is.vector(x) | !all(is.finite(x))) {
      stop("x should be a vector comprised of finite numbers")
   }
   # build det
   det <- det.construct(matrix(x, nrow = 1), mode = mode,
                        lb = bounds[1], ub = bounds[2], alphag = alpha,
                        progress = FALSE, cores = cores)
   leafs <- det.leafs(det)
   # plot det
   d <- matrix(NA, nrow = length(leafs$p), ncol = 2) # pdf estimate
   x <- d # x-axis
   for (k in 1:nrow(d)) { # loop over distribution elements or tree leafs
      x[k,] <- c(leafs$lb[k], leafs$lb[k]+leafs$size[k]) # element bounds
      # element distribution
      if (leafs$p[k] == 0) {d[k,] <- 0} else { # theta undefined for empty elements
         d[k,] <- leafs$p[k] * (c(-1,1)*leafs$theta[[k]]/2 + 1)
      }
   }
   m <- order(x[,1]) # sort elements for plotting
   x <- as.vector(t(x[m,])); d <- as.vector(t(d[m,]))
   graphics::plot(x, d, type = "l", ylim = c(0,max(d)), main = main)
}
