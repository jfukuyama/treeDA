#' Tree-based discriminant analysis
#'
#' A package for performing sparse, tree-based discriminant analysis.
#'
#' This package contains functions for building sparse,
#' tree-structured models for classification. The method is based on
#' the idea that when our predictors are structured according to a
#' tree, we can create an expanded feature space containing both the
#' original leaf predictors as well as node predictors, which
#' correspond to sums or averages across the leaves descending from
#' them. Without some sort of regularization this problem would be
#' unidentifiable, but with the regularization provided by sparse
#' discriminant analysis we get stable solutions.
#'
#' The package fits a sparse discriminant model in the expanded
#' feature space and translates the results back to the leaf space, so
#' that the interpretation can be purely in terms of the original
#' predictors. The package also includes functions to perform cross
#' validation to pick the sparsity level and plotting commands to
#' visualize the tree and the fitted coefficient vectors.
#'
#' The main function in this package is \code{\link{treeda}}, which
#' fits a sparse tree-based discriminant model. Additional functions
#' provided are \code{\link{treedacv}}, which performs
#' cross-validation to determine the correct sparsity level, and
#' functions to plot the resulting coefficient vectors along the tree
#' (\code{\link{plot_coefficients}}).
"_PACKAGE"
