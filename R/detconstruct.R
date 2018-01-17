#' Distribution Element Tree (DET) Construction
#'
#' The function \code{det.construct} generates a distribution element tree DET from available data. The DET can be used firstly in connection with \code{\link{det.query}} for density estimation. Secondly, with \code{\link{det.rnd}}, DETs can be used for smooth bootstrapping or more specifically conditional or unconditional random number generation.
#'
#' @param dta matrix with \code{d} rows representing components or dimensions and \code{n} columns corresponding to data points or samples.
#' @param mode order of distribution elements applied, default is \code{mode = 2}. Use \code{+/-1} for constant or \code{+/-2} for linear elements. \code{mode > 0} and \code{mode < 0} lead to equal-size and -score splits, respectively, in the element-refinement process.
#' @param lb,ub vectors of length \code{d} with lower and upper sample-space bounds. If not provided or set to \code{NA} or \code{0}, the bounds are determined from the data \code{dta}. If bounds are provided or given as \code{0}, the data is not pre-whitened before the DET is computed.
#' @param alphag,alphad significance levels for goodness-of-fit and independence tests, respectively, in element refinement or splitting process. Default is \code{alphag = alphad = 1.0e-3}. \code{alphad} is irrelevant for univariate data \code{dta} with \code{d = 1}.
#' @param progress optional logical, if set to \code{TRUE}, a progress report about the DET construction process is provided.
#' @param cores for large datasets, \code{cores > 1} allows for accelerated tree construction with parallel branch splitting. The default is \code{cores = 0}, which allocates half of the available cores (see \code{\link[parallel]{detectCores}}). \code{cores = 1} corresponds to serial tree construction.
#'
#' @return A DET object which describes the DET structure and pre-white transform is returned.
#' @references Meyer, D.W. (2016) \url{http://arxiv.org/abs/1610.00345} or Meyer, D.W., Statistics and Computing (2017) \url{https://doi.org/10.1007/s11222-017-9751-9} and Meyer, D.W. (2017) \url{http://arxiv.org/abs/1711.04632}
#' @export
#'
#' @examples
#' ## Gaussian data
#' require(stats)
#' det <- det.construct(t(rnorm(1e5)), cores = 2) # default linear elements
#' x <- t(seq(-4,4,0.05)); p <- det.query(det, x, cores = 2); plot(x, p, type = "l")
#'
#' ## piecewise uniform data with peaks
#' x <- matrix(c(rep(0,1e3),rep(1,1e3),2*runif(1e4),
#'               rep(0,5e2),rep(1,25e2),2*runif(9e3)), nrow = 2, byrow = TRUE)
#' det <- det.construct(x, mode = 1, lb = 0, ub = 0, cores = 2) # constant elements
det.construct <- function(dta, mode = 2, lb = NA, ub = NA, alphag = 1.0e-3, alphad = 1.0e-3,
                          progress = TRUE, cores = 0) {
   # setup cluster for parallel processing
   if (cores <= 0) {cores <- max(1,ceiling(parallel::detectCores()/2))}
   clster <- parallel::makeCluster(cores)
   # helper function to erase previously written progress report
   erasereport <- function() {cat(paste(rep("\b",53), collapse = ""))}
   # check data
   if (!is.matrix(dta)) {stop("invalid data, expecting d x n matrix with d dimensions and n samples")}
   if ((!all(is.finite(dta))) | is.null(dta) | (dim(dta)[2] == 0)) {
      stop("invalid data (NULL, NA, NaN or Inf), fix invalid entries")
   }
   d <- dim(dta)[1] # number of dimensions
   n <- dim(dta)[2] # number of samples
   # probability space or root cuboid bounds and principal axis transform
   if (is.na(lb[1])) { # initialize principal axes or pre-whitening transform
      if (dim(unique(dta, MARGIN = 2))[2] > 1) { # cov-matrix defined
         A <- t(eigen(stats::cov(t(dta)))$vectors); mu <- rowMeans(dta)
         dta <- A%*%(dta-mu%*%t(rep(1,n))) # transform to principal axes or pre-whitening
      } else { # cov-matrix undefined
         A <- diag(nrow = d); mu <- rep(0,d) # no pre-whitening
      }
      lb <- 0; ub <- 0 # trigger initialization of bounds
   } else { # lb exists -> no principal axes transform or pre-whitening
      A <- diag(nrow = d); mu <- rep(0,d)
   }
   if (all(lb == 0) & all(ub == 0)) { # initialize bounds if not given
      lb <- apply(dta,1,min); ub <- apply(dta,1,max)
      # treat singular case with lb = ub
      for (k in 1:d) {
         if (ub[k] == lb[k]) {
            ub[k] <- lb[k] + (abs(lb[k])+(lb[k] == 0)) * 10*.Machine$double.eps
         }
      }
   } else { # check if data is within bounds given
      if (!(all(lb%*%t(rep(1,n)) <= dta) & all(dta <= ub%*%t(rep(1,n))))) {
         stop("data is outside of bounds provided (lb <= dta <= ub)")
      }
   }
   dta <- (dta-lb%*%t(rep(1,n))) / ((ub-lb)%*%t(rep(1,n))) # normalize
   # initialize det fields (de[k] is given by tree[k,], sd[k], sp[k], p[k], & theta[[k]])
   mem <- ceiling(log(n)) # size of de memory block
   tree <- matrix(rep(NA,3*mem), ncol = 3); # tree[k,] contains indices of [parent child1 child2] of de[k]
   sd <- rep(NA,mem); sp <- rep(NA,mem) # split dimension and position of de leading to children
   p <- rep(NA,mem) # de probability densities
   theta <- rep(list(NA),mem) # de parameters
   # initialize lists with interim de data based on root de
   dei <- c(1) # list of interim elements
   sze <- list(rep(1,d)) # size of root de (normalized)
   x <- list(dta) # samples belonging to interim elements
   dt <- dimstosplit(x[[1]],sze[[1]],mode,alphag,alphad) # split of root de
   thetai <- list(dt$theta) # list of pdf parameters of interim elements
   tosplit <- list(dt$dim) # dimension(s) of interim elements to split (= NA for no split)

   # construct de tree level-by-level
   m <- 0 # number of tree leafs or final elements
   ne <- 1 # number of interim and final elements
   nt <- 0 # number of tree levels or tree depth
   vol <- 0 # volume of final elements
   while (length(dei) > 0) { # loop over tree levels

      # interim elements with tosplit = NA become final elements
      ind <- seq(along = tosplit)[is.na(tosplit)]
      for (k in ind) {vol <- vol + prod(sze[[k]])}
      # keep track of probability densities of final elements
      if (dim(tree)[1] < ne+2*sum(!is.na(tosplit))) { # allocate memory as needed
         k <- max(mem, ne+2*sum(!is.na(tosplit)) - dim(tree)[1]) # (required) memory for next level
         tree <- rbind(tree, matrix(rep(NA,3*k), ncol = 3))
         sd <- c(sd, rep(NA,k)); sp <- c(sp, rep(NA,k))
         p <- c(p, rep(NA,k)); theta <- c(theta, rep(list(NA),k))
      }
      for (k in ind) { # set probability densities of final elements
         p[dei[k]] <- (dim(x[[k]])[2])/n * 1/prod(sze[[k]]) # probability density of element
         theta[[dei[k]]] <- thetai[[k]] # pdf parameters of element
         m <- m+1 # counter of final elements
      }
      # discard final elements from list of interim elements
      ind <- seq(along = tosplit)[!is.na(tosplit)] # indices of interim elements
      dei <- dei[ind]; sze <- sze[ind]; x <- x[ind]
      tosplit <- tosplit[ind]; thetai <- thetai[ind]

      # get ready to split interim elements
      if (progress) {
         if (nt > 0) {erasereport()}
         cat(sprintf("%6.2f%% completed, %6i elements, tree depth %6i", 100*vol, m, nt))
      }
      if (length(dei) == 0) {next} # no interim elements left -> done
      nt <- nt + 1 # add new tree level
      # package data for parallel splitting of interim elements
      parents <- list(ind = 1:length(dei), tosplit = tosplit, x = x, sze = sze)
      parents <- lapply(seq_along(parents[[1]]), function(k) lapply(parents, "[[", k)) # transpose list
      # split interim elements in parallel (parent -> child1 + child2)
      cl <- parallel::parLapply(clster, parents, function(parent) {
         # split
         children <- de.split(parent$tosplit[1], parent$x, parent$sze, ne+2*parent$ind-1, mode)
         # prepare next splits
         if (length(parent$tosplit) > 1) { # 2nd split dimension from previous iteration
            children <- c(children, list(dim1 = parent$tosplit[2], theta1 = NA,
                                         dim2 = parent$tosplit[2], theta2 = NA))
         } else { # no 2nd split dimension -> de hypothesis testing
            dt <- dimstosplit(children$x1, children$sze1, mode, alphag, alphad)
            children <- c(children, list(dim1 = dt$dim, theta1 = dt$theta))
            dt <- dimstosplit(children$x2, children$sze2, mode, alphag, alphad)
            children <- c(children, list(dim2 = dt$dim, theta2 = dt$theta))
         }
         return(children)
      })
      # unpack data of children
      dei1 <- rep(NA,length(dei)); thetai1 <- rep(list(NA),length(dei)); tosplit1 <- thetai1; sze1 <- thetai1; x1 <- thetai1
      dei2 <- dei1; thetai2 <- thetai1; tosplit2 <- tosplit1; sze2 <- sze1; x2 <- x1
      pos <- rep(NA,length(dei))
      for (k in 1:length(dei)) {
         dei1[k] <- cl[[k]]$dei1; sze1[[k]] <- cl[[k]]$sze1; x1[[k]] <- cl[[k]]$x1
         dei2[k] <- cl[[k]]$dei2; sze2[[k]] <- cl[[k]]$sze2; x2[[k]] <- cl[[k]]$x2
         pos[k] <- cl[[k]]$pos
         tosplit1[[k]] <- cl[[k]]$dim1; thetai1[[k]] <- cl[[k]]$theta1
         tosplit2[[k]] <- cl[[k]]$dim2; thetai2[[k]] <- cl[[k]]$theta2
      }

      # update tree based on new elements (children)
      for (k in 1:length(dei)) {
         # parent
         ind <- dei[k]
         sd[ind] <- tosplit[[k]][1]
         sp[ind] <- pos[k]
         tree[ind,2:3] <- c(dei1[k],dei2[k]) # reference to children
         # children
         tree[dei1[k],1] <- ind; tree[dei2[k],1] <- ind # references to parent
      }
      ne <- ne+2*length(dei) # number of interim and final elements

      # prepare interim elements data for next tree level
      dei <- c(dei1,dei2); thetai <- c(thetai1,thetai2)
      tosplit <- c(tosplit1,tosplit2); sze <- c(sze1,sze2); x <- c(x1,x2)
   } # loop over tree levels

   if (progress) {erasereport()}
   parallel::stopCluster(clster)
   # discard unnecessary memory
   tree <- tree[1:ne,,drop = FALSE]; sd <- sd[1:ne]; sp <- sp[1:ne]; p <- p[1:ne]; theta <- theta[1:ne]
   # assemble output
   return(list(tree = tree, sd = sd, sp = sp, p = p, theta = theta, m = m, nt = nt,
               A = A, mu = mu, lb = lb, ub = ub))
}

#' Determine Split Dimension(s)
#'
#' Determine the split dimensions of an existing distribution element in the DET refinement-process based on statistical tests.
#'
#' @param x data in element given as matrix with \code{d} rows or components and \code{n} columns or samples.
#' @param sze vector representing the element size.
#' @param mode element order and split mode as detailed in \code{\link{det.construct}}.
#' @param alphag,alphad significance levels for goodness-of-fit and independence tests, respectively, in element refinement or splitting process. \code{alphad} is irrelevant for univariate data \code{x} with \code{d = 1}.
#'
#' @return An object comprised of the split dimension(s) or \code{NA} for no split, and the resulting child element parameters is returned.
dimstosplit <- function(x, sze, mode, alphag, alphad) {
   # is there any data?
   if (dim(x)[2] == 0) {return(list(dim = NA, theta = NA))}
   # proceed with data
   d <- dim(x)[1] # number of dimensions
   theta <- rep(0,d); p1 <- matrix(rep(0,2*d), nrow = d); h <- rep(FALSE,d)
   p1[,2] <- 1/sze # keep track of de side lengths
   for (k in 1:d) { # loop over dimensions
      # don't test based on just one unique sample, accept hypothesized distribution instead
      if (length(unique(x[k,])) <= 1) {next}
      # test distribution
      tst <- switch(abs(mode),
         c(chi2testuniform(x[k,], alphag), s = 0), # 1, hypothesis is uniform marginal distribution
         chi2testlinear(x[k,], alphag) # 2, hypothesis is linear marginal distribution
      )
      h[k] <- tst$h; p1[k,1] <- tst$p; theta[k] <- tst$s
   }
   dimens <- NA
   # was any marginal pdf test positive (null hypothesis rejected)?
   if (any(h)) { # yes, determine biggest deviating dimension (smallest p-value) and return
      # in case of equal p-value, sort by 1/size (equilibrate side lengths)
      m <- order(p1[,1], p1[,2])
      # select rejecting dimension
      for (k in 1:d) {
         if (h[m[k]]) {dimens <- m[k]; break()}
      }
   } else {
      if ((d > 1) & (dim(unique(x, MARGIN = 2))[2] >= 1)) { # no, check independence hypothesis
         tst <- chi2indeptest(x, alphad); h <- tst$h; p2 <- tst$p
         pl <- matrix(rep(0,d^2 - d), ncol = 2) # list of p-values
         dl <- pl # list of dependent dimensions
         nc <- 0 # counter of dependent dimension pairs
         for (j in 1:(d-1)) {
            for (k in (1+j):d) {
               if (h[j,k]) { # independence hypothesis is rejected
                  nc <- nc+1
                  pl[nc,1] <- p2[j,k] # p-value
                  pl[nc,2] <- 1/(sze[j]*sze[k]) # de side lengths
                  dl[nc,] <- c(j,k) # dimensions
               }
            }
         }
         # check if data is independent
         if (nc > 0) {
            # most dependent dimension at the top (smallest p-value)
            # in case of equal p-value, sort by 1/area (equilibrate side lengths)
            dimens <- order(pl[1:nc,1], pl[1:nc,2])
            dimens <- dl[dimens[1],]
            # smaller p-value first
            if (p1[dimens[1],1] >= p1[dimens[2],1]) {dimens <- rev(dimens)}
         }
      }
   }
   return(list(dim = dimens, theta = theta))
}

#' Split a Distribution Element
#'
#' Splits a parent distribution element characterized by \code{x}, \code{sze}, and \code{id} along dimension \code{dimens} based on \code{mode} into two child elements.
#'
#' @param dimens split dimension.
#' @param x data in element given as matrix with \code{d} rows or components and \code{n} columns or samples.
#' @param sze vector representing the size of the element.
#' @param id element index or identifier.
#' @param mode for splitting: \code{> 0} or \code{< 0} for equal-size or -score split, respectively.
#'
#' @return Object containing the properties of the resulting two child distribution elements and the split position within the parent element.
de.split <- function(dimens, x, sze, id, mode) {
   # determine position of subdivision
   if ((length(unique(x[dimens,])) <= 1) | (mode > 0)) { # empty or only one unique point
      pos <- 0.5
   } else { # multiple points
      pos <- stats::median(unique(x[dimens,]))
   }
   # first subcuboid
   de1 <- id
   sze1 <- sze; sze1[dimens] <- sze1[dimens]*pos
   # second subcuboid
   de2 <- id+1
   sze2 <- sze; sze2[dimens] <- sze2[dimens]*(1-pos)
   # samples in subcuboids
   if (dim(x)[2] == 0) { # no data
      x1 <- matrix(1, nrow = dim(x)[1], ncol = 0); x2 <- x1
   } else { # split data
      ind <- (x[dimens,] <= pos) # samples in first subcuboid
      x1 <- x[,ind,drop = FALSE]
      if (dim(x1)[2] > 0) {x1[dimens,] <- x1[dimens,]/pos} # rescale dim.
      x2 <- x[,!ind,drop = FALSE]
      if (dim(x2)[2] > 0) {x2[dimens,] <- (x2[dimens,]-pos)/(1-pos)} # rescale dim.
   }
   return(list(dei1 = de1, sze1 = sze1, x1 = x1, dei2 = de2, sze2 = sze2, x2 = x2, pos = pos))
}
