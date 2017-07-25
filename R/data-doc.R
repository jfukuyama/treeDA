#' Example dataset
#'
#' A small example dataset with three components, stored as a list
#' with a vector containing the classes (\code{response}), a matrix
#' containing the predictor variables (\code{predictors}), and a tree
#' describing the relationships between the predictor variables
#' (\code{tree}). The dataset consists of 50 samples divided into two
#' classes and 100 taxa/predictor variables, related to each other by
#' a random tree (generated with \code{ape::rtree}). A set of 42 taxa
#' descending from one internal node are all over-represented in one
#' class and under-represented in the other. The \code{predictors}
#' element in the list contains real numbers, not counts, and is
#' supposed to reflect normalized taxon abundances (e.g.,
#' normalization using the variance-stabilizing transformation in
#' DESeq2).
#'
#' @format A list containing response variables, predictor variables,
#'     and a tree describing the relationship between the predictor
#'     variables.
#' @name treeda_example
NULL
