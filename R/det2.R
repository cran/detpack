#' Density Estimation for Bivariate Data Based on Distribution Element Trees
#'
#' Two-dimensional piecewise linear or constant probability density estimator based on distribution element trees (DETs).
#'
#' @param x matrix with two rows containing data points (columns).
#' @param mode order of distribution elements applied, default is \code{mode = 2}. Use \code{+/-1} for constant or \code{+/-2} for linear elements. \code{mode > 0} and \code{mode < 0} lead to equal-size and -score splits, respectively, in the element-refinement process.
#' @param bounds \code{list(lb,ub)}, where \code{lb} and \code{ub} are vectors representing the lower and upper bounds of the probability space. If both bounds are set to NA (default) or 0, the bounds are determined based on the data \code{x}. Additionally, if the bounds are set to 0, pre-whitening is not applied to the data.
#' @param alphag,alphad significance levels for goodness-of-fit and independence tests, respectively, in element refinement or splitting process. Default is \code{alphag = alphad = 1.0e-3}.
#' @param main an overall plot title, see \code{\link[graphics]{title}}.
#' @param np vector with number of pixels in the first and second probability-space directions used for rendering the DET-based probability density estimate (default is \code{100x100}).
#' @param cores number of cores for parallel tree construction and query (rendering). Default is \code{0} for automatic allocation, see \code{cores} in \code{\link{det.construct}} and \code{\link{det.query}}.
#'
#' @export
#'
#' @examples
#' ## uniform
#' require(stats)
#' det2(rbind(runif(1e4),1+2*runif(1e4)), mode = 1, bounds = list(c(-0.1,0),c(1.1,4)), cores = 1)
#' det2(rbind(1:100,101:200), mode = 2, np = c(160,160), cores = 1) # data on a line
#'
#' ## Gaussian
#' require(stats); require(graphics)
#' n <- 1e3; x <- rnorm(n)
#' x <- matrix(c(x, x+rnorm(n,0,0.5)), ncol = 2)
#' split.screen(c(2,2))
#' screen(3); plot(x, type = "p", pch = ".", main = "data")
#' screen(1); det2(t(x), mode = 1, cores = 1, main = "constant det estimator")
#' screen(2); det2(t(x), cores = 1, main = "linear det estimator")
#' screen(4); det2(t(x), mode = 1, bounds = list(0,0), cores = 2, main = "const. det, no pre-white")
det2 <- function(x, mode = 2, bounds = list(NA,NA), alphag = 1.0e-3, alphad = 1.0e-3,
                 main = NULL, np = c(100,100), cores = 0) {
   # check data
   errorstrg <- "x should be a matrix with two rows comprised of finite numbers"
   if (is.matrix(x) & all(is.finite(x))) {
      if (dim(x)[1] != 2) {stop(errorstrg)}
   } else {stop(errorstrg)}
   # build det
   det <- det.construct(x, mode = mode, lb = bounds[[1]], ub = bounds[[2]],
                        alphag = alphag, alphad = alphad,
                        progress = FALSE, cores = cores)
   # bounds
   if (is.na(bounds[1])) {bounds[[1]] = apply(x, 1, min)}
   else {if (all(bounds[[1]] == 0)) {bounds[[1]] = det$lb}}
   if (is.na(bounds[2])) {bounds[[2]] = apply(x, 1, max)}
   else {if (all(bounds[[2]] == 0)) {bounds[[2]] = det$ub}}
   # plot or render det using grid points xp
   x1 <- bounds[[1]][1] + (bounds[[2]][1]-bounds[[1]][1])*((1:np[1])-0.5)/np[1]
   x2 <- bounds[[1]][2] + (bounds[[2]][2]-bounds[[1]][2])*((1:np[2])-0.5)/np[2]
   xp <- rbind(rep(x1,np[2]), rep(x2,each = np[1]))
   y <- det.query(det, xp, cores = cores); y <- matrix(y, nrow = np[1])
   graphics::image(x = x1, y = x2, z = y, col = grDevices::gray((100:0)/100), main = main)
}
