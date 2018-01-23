#' Goodness-of-Fit Test for Uniform Distribution
#'
#' Pearson's Chi-square test for goodness-of-fit of a uniform distribution.
#'
#' @param x vector with normalized data between 0 and 1.
#' @param alpha significance level.
#'
#' @return Object with test outcome \code{h = TRUE/FALSE} meaning rejection/acceptance of uniform null hypothesis, and p-value or confidence level of acceptance.
#'
# examples
# require(stats)
# n <- 1e5; lb <- 4; ub <- 4.01
# ## uniformly distributed data
# x <- lb + (ub-lb)*runif(n) # data in [0,1]
# tst <- chi2testuniform((x-lb)/(ub-lb), 1.0e-3)
# ## clipped Gaussian data
# x <- rnorm(n,(lb+ub)/2,0.1*(ub-lb)); x[x<lb | x>ub] <- (lb+ub)/2
# tst <- chi2testuniform((x-lb)/(ub-lb), 1.0e-3)
chi2testuniform <- function(x, alpha) {
   n <- length(x) # number of samples
   # binning
   tbl <- chi2testtable(x, alpha) # observed count and bin edges
   k <- length(tbl$co) # number of bins
   if (k <= 1) {return(list(h = FALSE, p = 0))} # no testing with just one bin
   ce <- n * (tbl$be[2:(k+1)]-tbl$be[1:k]) # expected count
   # perform composite Pearson's chi^2 test
   chi2 <- sum(((tbl$co-ce)^2)/ce) # test statistic
   p <- stats::pchisq(chi2, k-1, lower.tail = FALSE) # p-value = 1 - chi^2 cdf
   h <- (p <= alpha) # accept or reject (FALSE or TRUE)
   return(list(h = h, p = p))
}

#' Goodness-of-Fit Test for Linear Distributions
#'
#' Composite Pearson's Chi-square test for goodness-of-fit of linear distributions.
#'
#' @param x vector with normalized data between 0 and 1.
#' @param alpha significance level.
#'
#' @return Object with test outcome \code{h = TRUE/FALSE} meaning rejection/acceptance of linear null hypothesis, and p-value or confidence level of acceptance.
#'
# examples
# require(stats)
# n <- 1e5; lb <- 4; ub <- 4.01
# ## linearly distributed data
# x <- runif(n); y <- runif(n); x <- x[x/2+1/2 >= y]
# x <- lb + (ub-lb)*x # data in [0,1]
# tst <- chi2testlinear((x-lb)/(ub-lb), 1.0e-3)
# ## clipped Gaussian data
# x <- rnorm(n,(lb+ub)/2,0.1*(ub-lb))
# x[x<lb | x>ub] <- (lb+ub)/2 # clipping to [0,1]
# tst <- chi2testlinear((x-lb)/(ub-lb), 1.0e-3)
chi2testlinear <- function(x, alpha) {
   cdf <- function(x,slp) return((2+slp*(x-1))*x / 2) # cdf linear distrib.
   n <- length(x) # number of samples
   # estimate mean to determine slope parameter of linear distribution
   m <- mean(x)
   if (m == 1/2) {s <- 0} # in case of var(x) = 0
   else {
      s <- 6*(2*m-1)
      # correction factor cf to minimize mean square error
      cf <- 1; if (length(x) > 1) {cf <- (n*s^2)/(n*s^2 + 4*36*stats::var(x))}
      s <- cf*s
   }
   s <- min(2, max(-2,s)) # limit s to make sure pdf >= 0
   # binning
   tbl <- chi2testtable(x, alpha) # observed count and bin edges
   k <- length(tbl$co) # number of bins
   if (k <= 1) {return(list(h = FALSE, p = 0, s = 0))} # no testing with just one bin
   ce <- rep(0,k) # expected count
   for (j in 1:k) {ce[j] <- cdf(tbl$be[j+1],s)}
   ce[2:k] <- ce[2:k] - ce[1:(k-1)]; ce <- n*ce
   # perform composite Pearson's chi^2 test
   chi2 <- sum(((tbl$co-ce)^2)/ce) # test statistic
   p <- stats::pchisq(chi2, k-2, lower.tail = FALSE) # p-value = 1 - chi^2 cdf
   h <- (p <= alpha) # accept or reject (FALSE or TRUE)
   return(list(h = h, p = p, s = s))
}

#' Pairwise Mutual Independence Test
#'
#' Pearson's Chi-square test of pairwise mutual independence.
#'
#' @param x vector with normalized data between 0 and 1.
#' @param alpha significance level.
#'
#' @return Object with test outcome \code{h = TRUE/FALSE} meaning rejection/acceptance of pairwise independence null hypothesis, and p-value or confidence level of acceptance.
#'
# examples
# ## x[1,] and x[3,] are independent, but x[2,] = (x[1,] + uniform noise)/2
# n <- 100; x <- matrix(stats::runif(3*n),nrow=3,ncol=n); x[2,] <- (x[1,]+x[2,])/2
# tst = chi2indeptest(x, 1.0e-3)
chi2indeptest <- function(x, alpha) {
   # initialize things
   n <- nrow(x) # number of datasets in x
   p <- matrix(rep(0,n*n), nrow = n, ncol = n) # p-value matrix
   h <- matrix(rep(TRUE,n*n), nrow = n, ncol = n) # matrix with test results
   # test for mutual independence
   for (i in 1:(n-1)) {
      for (j in (i+1):n) { # upper triangle
         # marginal observed counts & bin edges
         tbl1 <- chi2testtable(x[i,], alpha, TRUE)
         tbl2 <- chi2testtable(x[j,], alpha, TRUE)
         n1 <- length(tbl1$co); n2 <- length(tbl2$co)
         # no need to test if there is just one bin
         if ((n1 == 1) | (n2 == 1)) {
            h[i,j] <- FALSE; h[j,i] <- FALSE # accept
            next
         }
         # contingency table with observed counts
         i1 <- cut(x[i,], tbl1$be, include.lowest = TRUE)
         i2 <- cut(x[j,], tbl2$be, include.lowest = TRUE)
         co <- table(i1,i2)
         # expected counts
         ce <- (tbl1$co %*% (t(tbl2$co)))/ncol(x)
         # chi^2 independence test
         chi2 <- sum(sum(((co-ce)^2)/ce)) # test statistic
         p[i,j] <- stats::pchisq(chi2, (n1-1)*(n2-1), lower.tail = FALSE) # p-value = 1 - chi^2 cdf
         p[j,i] <- p[i,j]
         h[i,j] <- (p[i,j] <= alpha); h[j,i] <- h[i,j] # accept or reject (FALSE or TRUE)
      }
   }
   return(list(h = h, p = p))
}

#' (Contingency) Tables for Chi-square Tests
#'
#' (Contingency) tables for Pearson's Chi-square goodness-of-fit and independence tests.
#'
#' @param x vector with normalized data between 0 and 1.
#' @param alpha significance level.
#' @param cf flag, if \code{TRUE}, bin counts for a multi-variate contingency table are applied (Bagnato et al., 2012). Otherwise rules for univariate tables are used (Cochran, 1952; Mann and Wald, 1942).
#' @references
#' Cochran, W.G., The Chi Square Test of Goodness of Fit. The Annals of Mathematical Statistics, 1952. 23(3): p. 315-345.
#'
#' Mann, H.B. and A. Wald, On the Choice of the Number of Class Intervals in the Application of the Chi Square Test. The Annals of Mathematical Statistics, 1942. 13(3): p. 306-317.
#'
#' Bagnato, L., A. Punzo, and O. Nicolis, The autodependogram: a graphical device to investigate serial dependences. Journal of Time Series Analysis, 2012. 33(2): p. 233-254.
#'
#' @return Object with observed counts \code{co} of data \code{x} inside bins defined by edges \code{be}.
#'
# examples
# x <- rnorm(100); x <- x-min(x); x <- x/max(x)
# tbl <- chi2testtable(x, 0.01)
chi2testtable <- function(x, alpha, cf = FALSE) {
   # first try without eliminating duplicates (unique is an expensive operation)
   xu <- x
   for (t in 1:2) {
      # number of bins based on rule of thumb [Cochran, Ann. of Math. Stat., 23(3),
      # section 7] and expression in [Mann & Wald, Ann. of Math. Stat., 13(3)], ...
      n <- length(xu) # number of (unique) samples
      c <- stats::qnorm(alpha, lower = FALSE)
      ks <- n/5
      kp <- 4*(2*(n-1)^2/c^2)^(1/5)
      # ... for contingency tables see [Bagnato, Punzo, & Nicolis, J. Time Ser. Anal., 33]
      if (cf) {ks <- sqrt(ks); kp <- sqrt(kp)}
      k <- max(1, floor(min(ks,kp))) # number of bins
      # put data into bins
      xu <- sort(xu)
      be <- c(rep(0, k),1) # bin edges
      for (j in 1:k) { # loop over bins
         ns <- floor(j/k*n+0.5) # index of last sample in current bin
         # floor is used instead of round because round(0.5) -> 0 not 1 (!= matlab)
         if (j < k) {be[j+1] <- (xu[ns+1] + xu[ns])/2} # all but last bins
      }
      # check if there are no zero bins
      if (all(be[2:length(be)] != be[1:(length(be)-1)])) {break}
      # second attempt, now without duplicates
      xu <- unique(x)
   }
   rm(xu) # free memory
   co <- cut(x, breaks = be, include.lowest = TRUE) # convert data to factor
   co <- table(co) # table with counts of factor levels
   return(list(co = co, be = be))
}
