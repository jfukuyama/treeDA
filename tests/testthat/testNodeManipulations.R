context("tree.manipulation")
library(ape)
library(treeDA)
test_that("Edges to children works", {
    tree = stree(n = 2, type = "balanced")
    etc = treeDA:::edgesToChildren(tree$edge)
    # node 3 is the root, its children are 1 and 2
    expect_equal(1 %in% etc[[3]], TRUE)
    expect_equal(2 %in% etc[[3]], TRUE)
    expect_equal(length(etc[[3]]), 2)
    # node 1 and 2 are leaf nodes and have no children
    expect_equal(length(etc[[1]]), 0)
    expect_equal(length(etc[[2]]), 0)
})

test_that("Predictors to nodes works", {
    tree = stree(n = 2, type = "balanced")
    leafPredictors = matrix(c(0,1,2,1), nrow = 2)
    fullPredictors = treeDA:::makeNodeAndLeafPredictors(leafPredictors, tree)
    # fullPredictors[,c(1,2)] should match leafPredictors
    expect_equal(sum(abs(fullPredictors[,c(1,2)] - leafPredictors)), 0)
    # the root node predictor should be equal to c(2,2)
    expect_equal(sum(abs(fullPredictors[,3] - c(2,2))), 0)
})


test_that("Make descendant matrix works", {
    tree = stree(n = 2, type = "balanced")
    dm = treeDA:::makeDescendantMatrix(tree)
    expect_equal(nrow(dm), 2)
    expect_equal(ncol(dm), 3)
    expect_equal(dm[1,1], 1)
    expect_equal(dm[2,2], 1)
    expect_equal(dm[1,2], 0)
    expect_equal(dm[2,1], 0)
    expect_equal(dm[1,3], 1)
    expect_equal(dm[2,3], 1)
})

test_that("Coefficient vector mapping works", {
    tree = stree(n = 2, type = "balanced")
    leafCoefficients = nodeToLeafCoefficients(c(0,1,2), tree)
    expect_equal(leafCoefficients[1], 2)
    expect_equal(leafCoefficients[2], 3)
})
