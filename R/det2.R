#' Density Estimation for Bivariate Data Based on Distribution Element Trees
#'
#' Two-dimensional piecewise linear or constant probability density estimator based on distribution element trees (DETs).
#'
#' @param dta matrix with two rows containing data (samples in columns).
#' @param mode order of distribution elements applied, default is \code{mode = 2}. Use \code{+/-1} for constant or \code{+/-2} for linear elements. \code{mode > 0} and \code{mode < 0} lead to equal-size and -score splits, respectively, in the element-refinement process.
#' @param bounds \code{list(lb,ub)}, where \code{lb} and \code{ub} are vectors representing the lower and upper bounds of the probability space. If both bounds are set to NA (default) or 0, the bounds are determined based on the data \code{dta}. Additionally, if the bounds are set to 0, pre-whitening is not applied to the data.
#' @param alphag,alphad significance levels for goodness-of-fit and independence tests, respectively, in element refinement or splitting process. Default is \code{alphag = alphad = 1.0e-3}.
#' @param main an overall plot title, see \code{\link[graphics]{title}}. If \code{main = NULL} (default), the density range is provided as a title.
#' @param nc number of contour levels (default is 20).
#' @param dtalim allows to limit the number of samples used in tests guiding the element splitting process. Default is \code{dtalim = Inf}, which corresponds to using all available samples, see \code{\link{det.construct}}.
#' @param cores number of cores for parallel tree construction. Default is \code{cores = 1} for serial processing, see \code{cores} in \code{\link{det.construct}}.
#'
#' @export
#'
#' @examples
#' ## uniform
#' require(stats)
#' det2(rbind(runif(5e3),1+2*runif(5e3)), mode = 1, bounds = list(c(-0.1,0),c(1.1,4)))
#' det2(rbind(1:100,101:200+runif(100)), mode = 2) # data on a line
#'
#' ## Gaussian
#' require(stats); require(graphics); require(grDevices)
#' n <- 5e3; x <- rnorm(n)
#' x <- matrix(c(x, x+rnorm(n,0,0.5)), ncol = 2)
#' split.screen(c(2,2))
#' screen(3); plot(x, type = "p", pch = ".", main = "data")
#' screen(1); det2(t(x), mode = 1, main = "constant det estimator")
#' screen(2); det2(t(x), main = "linear det estimator")
#' screen(4); det2(t(x), mode = 1, bounds = list(0,0), main = "const. det, no pre-white")
det2 <- function(dta, mode = 2, bounds = list(NA,NA), alphag = 1.0e-3, alphad = 1.0e-3,
                 main = NULL, nc = 20, dtalim = Inf, cores = 1) {
   # check data
   errorstrg <- "dta should be a matrix with two rows comprised of finite numbers"
   if (is.matrix(dta) & all(is.finite(dta))) {
      if (nrow(dta) != 2) {stop(errorstrg)}
   } else {stop(errorstrg)}
   # build det
   det <- det.construct(dta, mode = mode, lb = bounds[[1]], ub = bounds[[2]],
                        alphag = alphag, alphad = alphad,
                        progress = FALSE, dtalim = dtalim, cores = cores)
   # bounds
   if (is.na(bounds[1])) {bounds[[1]] = apply(dta, 1, min)}
   else {if (all(bounds[[1]] == 0)) {bounds[[1]] = det$lb}}
   if (is.na(bounds[2])) {bounds[[2]] = apply(dta, 1, max)}
   else {if (all(bounds[[2]] == 0)) {bounds[[2]] = det$ub}}
   # determine pdf range
   leafs <- det.leafs(det)
   d <- matrix(rep(0,length(leafs$p)*4), nrow = length(leafs$p)) # pdf element corners
   for (k in 1:length(leafs$p)) { # loop over distribution elements or tree leafs
      # element distribution
      if (leafs$p[k] == 0) { # theta undefined for empty elements
         d[k,] <- 0
      } else { # theta defined
         d[k,1] <- leafs$p[k] * (-leafs$theta[[k]][1]/2 + 1) * (-leafs$theta[[k]][2]/2 + 1)
         d[k,2] <- leafs$p[k] * (+leafs$theta[[k]][1]/2 + 1) * (-leafs$theta[[k]][2]/2 + 1)
         d[k,3] <- leafs$p[k] * (+leafs$theta[[k]][1]/2 + 1) * (+leafs$theta[[k]][2]/2 + 1)
         d[k,4] <- leafs$p[k] * (-leafs$theta[[k]][1]/2 + 1) * (+leafs$theta[[k]][2]/2 + 1)
      }
   }
   dmin <- 0; dmax <- max(d)
   # plot det element-by-element
   if (is.null(main)) {
      main <- paste("white",format(dmin, digits = 3),
                    "...", format(dmax, digits = 3, scientific = TRUE), "black")
   } # plot title = density range if main is not given
   graphics::plot(0, type="n", xlab="x1", ylab="x2", main = main,
                  xlim=c(bounds[[1]][1], bounds[[2]][1]),
                  ylim=c(bounds[[1]][2], bounds[[2]][2]))
   xy <- matrix(rep(0,2*4), nrow = 2) # element corner points
   for (k in 1:length(leafs$p)) { # loop over distribution elements or tree leafs
      xy[,1] <- leafs$lb[,k]
      xy[,2] <- c(leafs$lb[1,k]+leafs$size[1,k], leafs$lb[2,k])
      xy[,3] <- c(leafs$lb[1,k]+leafs$size[1,k], leafs$lb[2,k]+leafs$size[2,k])
      xy[,4] <- c(leafs$lb[1,k], leafs$lb[2,k]+leafs$size[2,k])
      xy <- t(det$A) %*% xy + det$mu %*% t(rep(1,4)) # pre-white transform
      contourRect(xy, d[k,], -abs(nc), dmin, dmax) # draw element
   }
}

#' Draw Contours in a Rectangle
#'
#' The function \code{contourRect} draws the z contour levels of a rectangular domain in x-y-space with z-values given at the corners of the rectangle.
#'
#' @param xy matrix with two rows and four columns containing x- and y-coordinates of the four corner points of the rectangle. The corner points are ordered in clockwise or counter-clockwise direction.
#' @param z vector with four z-values at the four corner points.
#' @param n \code{abs(n)} gives the number of local or global contour levels. If \code{n > 0}, \code{n} local contours are drawn within \code{[min(z),max(z)]}. If \code{n < 0}, \code{n} global contours result in \code{[zlb,zub]}, but only the contours falling inside \code{[min(z),max(z)]} are drawn.
#' @param zlb,zub determines the global range of z-values used to determine the contour colors. All values in \code{z} have to be contained in \code{[zlb,zub]}.
#'
# examples
# graphics::plot(0, type="n", xlab="", ylab="", xlim=c(0, 3), ylim=c(10, 13))
# zlb <- -0.5; zub <- 2; n <- -8
# contourRect(matrix(c(0,10,1,10,1,11.5,0,11.5), nrow = 2), c(2,1,1.6,1.3), n, zlb, zub)
# contourRect(matrix(c(0,11.5,1,11.5,1,13,0,13), nrow = 2), c(1.3,1.6,2,2), n, zlb, zub)
# contourRect(matrix(c(1,10,3,10,3,11.5,1,11.5), nrow = 2), c(1,0.2,2,1.6), n, zlb, zub)
# contourRect(matrix(c(1,11.5,3,11.5,3,13,1,13), nrow = 2), c(1.6,2,2,2), n, zlb, zub)
contourRect <- function(xy, z, n = 20, zlb = 0, zub = 1) {
   #color <- function(v) {return(grDevices::hsv(1,0,max(0,min(v,1))))} # gray scale
   color <- function(v) {
      v <- 1-v
      r <- stats::spline(x = 1:9, y = c(0,0.15,0.3,0.6,1,0.9,0.9,0.9,1), xout = 9*v)$y
      g <- stats::spline(x = 1:9, y = c(0,0.15,0.15,0.2,0.25,0.5,0.75,0.9,1), xout = 9*v)$y
      b <- stats::spline(x = 1:9, y = c(0,0.5,0.75,0.5,0.15,0,0.1,0.5,1), xout = 9*v)$y
      return(grDevices::rgb(max(0,min(r,1)),max(0,min(g,1)),max(0,min(b,1))))
   } # color scale
   zmin <- min(z); zmax <- max(z)
   # check input
   if ((zlb > zmin) | (zub < zmax)) stop("z should be inside [zlb,zub]")
   # in case of global contours, determine local zmin, zmax & n
   if (n < 0) {
      dz <- (zub-zlb)/(-n+1) # contour level increment
      zmin <- zlb + floor((zmin-zlb)/dz)*dz
      zmax <- zlb + ceiling((zmax-zlb)/dz)*dz
      n <- ceiling((zmax-zlb)/dz)-floor((zmin-zlb)/dz)-1 # number of contour levels to draw
      rm(dz)
   }
   # draw lowest contour level at zmin (rectangle)
   graphics::polygon(xy[1,], xy[2,], border = NA,
                     col = color((zmin-zlb)/(zub-zlb)))
   if ((zmax > zmin) & (n != 0)) {
      # close rectangle
      xy <- cbind(xy,xy[,1]); z <- c(z, z[1])
      # loop over levels (incl. zmax)
      for (l in 1:(n+1)) {
         zl <- l*(zmax-zmin)/(n+1) + zmin # level l
         # determine m polygon points of level l
         xyp <- matrix(rep(0,2*6), nrow = 2); m <- 0 # polygon points
         zc <- (c(z, z[1]) >= zl) # z of corner points in polygon
         for (p in 1:4) { # loop over edges connecting corners p & p+1
            if (zc[p]) {m <- m+1; xyp[,m] <- xy[,p]}
            if (xor(zc[p],zc[p+1])) {
               m <- m+1
               xyp[,m] <- xy[,p] +
                  (xy[,p+1]-xy[,p]) * (zl-z[p])/(z[p+1]-z[p])
            }
         } # for p
         # draw contour polygon (involving more than two points)
         if (m > 2) {
            graphics::polygon(xyp[1,1:m], xyp[2,1:m], border = NA,
                              col = color((zl-zlb)/(zub-zlb)))
         }
      } # for l
   }
}

