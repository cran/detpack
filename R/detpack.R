#' Distribution Element Trees for Density Estimation and Bootstrapping
#'
#' Distribution element trees (DETs) enable the estimation of probability densities based on (possibly large) datasets.
#' Moreover, DETs can be used for random number generation or smooth bootstrapping both in unconditional and conditional modes.
#' In the latter mode, information about certain probability-space components is taken into account when sampling the remaining probability-space components.
#'
#' The function \code{\link{det.construct}} translates a dataset into a DET.
#' To evaluate the probabilty density based on a precomputed DET at arbitrary query points, \code{\link{det.query}} is used.
#' The specific functions \code{\link{det1}} and \code{\link{det2}} provide density evaluation and plotting for one- and two-dimensional cases.
#' (Un)conditional smooth bootstrapping from an available DET, can be performed by \code{\link{det.rnd}}.
#' To inspect the structure of a DET, the functions \code{\link{det.de}} and \code{\link{det.leafs}} are useful. While \code{\link{det.de}} enables the extraction of an individual distribution element from the tree, \code{\link{det.leafs}} extracts all leaf elements at branch ends.
#'
#' @references
#' Distribution element tree basics and density estimation, see Meyer, D.W. (2016) \url{http://arxiv.org/abs/1610.00345} or Meyer, D.W., Statistics and Computing (2017) \url{https://doi.org/10.1007/s11222-017-9751-9}.
#'
#' DETs for smooth bootstrapping, see Meyer, D.W. (2017) \url{http://arxiv.org/abs/1711.04632}.
#'
#' @author Daniel Meyer, \email{meyerda@@ethz.ch}
#'
#' @docType package
#' @name detpack
NULL
