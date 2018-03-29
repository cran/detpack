#' Density Estimation for Univariate Data Based on Distribution Element Trees
#'
#' One-dimensional piecewise linear or constant probability density estimator based on distribution element trees (DETs).
#'
#' @param dta vector with data
#' @param mode order of distribution elements applied, default is \code{mode = 2}. Use \code{+/-1} for constant or \code{+/-2} for linear elements. \code{mode > 0} and \code{mode < 0} lead to equal-size and -score splits, respectively, in the element-refinement process.
#' @param bounds \code{c(lb,ub)}, where \code{lb} and \code{ub} are lower and upper bounds of the probability space. If both bounds are set to 0 (default), the bounds are determined based on the data \code{dta}.
#' @param alpha significance level for goodness-of-fit testing in element refinement or splitting process. Default is \code{alpha = 1.0e-3}.
#' @param main an overall plot title, see \code{\link[graphics]{title}}.
#' @param dtalim allows to limit the number of samples used in tests guiding the element splitting process. Default is \code{dtalim = Inf}, which corresponds to using all available samples, see \code{\link{det.construct}}.
#' @param cores number of cores for parallel tree construction. Default is \code{1} for serial construction, see \code{\link{det.construct}}.
#'
#' @export
#'
#' @examples
#' require(stats)
#' det1(rbeta(5e5, shape1 = 1.05, shape2 = 0.8), mode = -1,
#'      bounds = c(-0.1,1.1), main = "beta, const. elements, equal-scores splits")
#' x <- seq(-0.1,1.1,0.005); lines(x, dbeta(x,shape1 = 1.05,shape2 = 0.8), col = "red")
#' det1(rbeta(5e5, shape1 = 1.05, shape2 = 0.8), mode = -2,
#'      bounds = c(-0.1,1.1), main = "beta, linear elements, equal-scores splits")
#' x <- seq(-0.1,1.1,0.005); lines(x, dbeta(x,shape1 = 1.05,shape2 = 0.8), col = "red")
#' det1(rnorm(5e5), mode = 2, cores = 1, main = "Gaussian, linear elements, equal-size splits")
#' x <- seq(-5,5,0.05); lines(x, dnorm(x), col = "red")
#' det1(runif(5e5), mode = 1, bounds = c(-0.1,1.1),
#'      main = "uniform, const. elements, equal-size splits")
#' x <- seq(-0.1,1.1,0.005); lines(x, dunif(x), col = "red")
det1 <- function(dta, mode = 2, bounds = c(0,0), alpha = 1.0e-3, main = NULL, dtalim = Inf, cores = 1) {
   # check data
   if (!is.vector(dta) | !all(is.finite(dta))) {
      stop("dta should be a vector comprised of finite numbers")
   }
   # build det
   det <- det.construct(matrix(dta, nrow = 1), mode = mode,
                        lb = bounds[1], ub = bounds[2], alphag = alpha,
                        progress = FALSE, dtalim = dtalim, cores = cores)
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
