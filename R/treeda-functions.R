#' Check predictors
#'
#' Checks whether the predictors are consistent with the tree
#' structure.
#' @keywords internal
checkPredictorsAndTree <- function(predictors, tree) {
    if(ncol(predictors) != length(tree$tip.label)) {
        stop("Predictor matrix has different number of columns than tree")
    }
    if(!all(colnames(predictors) == tree$tip.label)) {
        stop("Column names of predictor matrix do not match the tip labels in the\ntree, you should modify the predictor matrix so that the columns\nare in the same order as the leaves of the tree and have matching\nnames")
    }
    if(is.null(colnames(predictors))) {
        warning("Predictor matrix has no column names, you should ensure that the\norder of the columns is the same as the order of the leaves of the\ntree")
    }
    return(0)
}

#' Makes a hash table with nodes and their children
#'
#' Takes the edge matrix from a phylo-class object and turns it into a
#' list where the entries are nodes and the elements are vectors with
#' the children of the nodes.
#' @keywords internal
edgesToChildren <- function(edges) {
    ## ape requires that the tips are numbered 1:n, nodes are numbered
    ## n+1:n+m (n tips, m internal nodes), root is numbered n+1. The
    ## first column of the edge matrix contains only the internal nodes
    etc = vector("list", max(edges))
    for(i in 1:nrow(edges)) {
        parent = edges[i,1]
        child = edges[i,2]
        if(is.null(etc[[parent]]))
            etc[[parent]] = child
        else
            etc[[parent]] = c(etc[[parent]], child)
    }
    return(etc)
}

#' Tree-based sparse discriminant analysis
#' 
#' Performs tree-structured sparse discriminant analysis using an
#' augmented predictor matrix with additional predictors corresponding
#' to the nodes and then translating the parameters back in terms of
#' only the leaves.
#'
#' @param response A factor or character vector giving the class to be
#'     predicted.
#' @param predictors A matrix of predictor variables corresponding to
#'     the leaves of the tree and in the same order as the leaves of
#'     the tree.
#' @param tree A tree of class \code{phylo}.
#' @param p The number of predictors to use.
#' @param k The number of components to use.
#' @param center Center the predictor variables?
#' @param scale Scale the predictor variables?
#' @param class.names Optional argument giving the class names.
#' @param check.consist Check consistency of the predictor matrix and
#'     the tree.
#' @param A A matrix describing the tree structure. If it has been
#'     computed before it can be passed in here and will not be
#'     recomputed.
#' @param ... Additional arguments to be passed to sda
#' @return An object of class \code{treeda}. Contains the coefficients
#'     in the original predictor space (\code{leafCoefficients}), the
#'     number of predictors used in the node + leaf space
#'     (\code{nPredictors}), number of leaf predictors used
#'     (\code{nLeafPredictors}), the projections of the samples onto
#'     the discriminating axes (\code{projections}), and the sparse
#'     discriminant analysis object that was used in the fit
#'     (\code{sda}).
#' @importFrom sparseLDA sda
#' @importFrom stats sd
#' @importFrom Matrix colMeans
#' @examples
#' data(treeda_example)
#' out.treeda = treeda(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     p = 1)
#' out.treeda
#' @export
treeda <- function(response, predictors, tree, p, k = nclasses - 1,
                   center = TRUE, scale = TRUE, class.names = NULL, check.consist = TRUE,
                   A = NULL, ...) {
    out = list()
    class(out) = "treeda"
    nclasses = length(unique(response))
    if(check.consist) {
        checkPredictorsAndTree(predictors, tree)
    }
    if(is.null(A))
        A = makeDescendantMatrix(tree)
    ep = as(predictors %*% A, "matrix")
    responseMatrix = makeResponseMatrix(response, class.names)
    if(center) {
        out$means = colMeans(ep)
        ep = sweep(ep, 2, STATS = out$means, FUN = "-")
    }
    if(scale) {
        out$sds = apply(ep, 2, sd)
        out$sds[out$sds == 0] = 1
        ep = sweep(ep, 2, STATS = out$sds, FUN = "/")
    }
    sda.out = sparseLDA::sda(x = ep, y = responseMatrix, stop = -p, ...)
    leafCoefficients = makeLeafCoefficients(sda.out, A, out$means, out$sds)
    out$leafCoefficients = leafCoefficients
    out$input = list(response = response, predictors = predictors, tree = tree)
    out$nPredictors = length(sda.out$varIndex)
    out$nLeafPredictors = sum(apply(out$leafCoefficients$beta, 1, function(x) any(x != 0)))
    out$sda = sda.out
    out$class.names = colnames(responseMatrix)

    out$projections = ep[,sda.out$varIndex] %*% sda.out$beta
    out$classProperties = makeClassProperties(response, out$projections)
    predictions = predict.treeda(out, newdata = predictors, newresponse = response,
        check.consist = FALSE)
    out$predictedClasses = predictions$classes
    out$predictionError = predictions$predictionError
    out$rss = predictions$rss
    return(out)
}

#' treeda cross validation
#'
#' Performs cross-validation of a \code{\link{treeda}} fit.
#'
#' @param response The classes to be predicted. 
#' @param predictors A matrix of predictors corresponding to the tips
#' of the tree.
#' @param tree A tree object of class \code{phylo}.
#' @param folds Either a single number corresponding to the number of
#' folds of cross-validation to perform or a vector of integers
#' ranging from 1 to the number of folds desired giving the partition
#' of the dataset.
#' @param pvec The values of p to use.
#' @param k The number of discriminating axes to keep.
#' @param center Center the predictors?
#' @param scale Scale the predictors?
#' @param class.names A vector giving the names of the classes.
#' @param ... Additional arguments to be passed to \code{\link{treeda}}. 
#'
#' @return A list with the value of p with minimum cv error
#'     (\code{p.min}), the minimum value of p with in 1 se of the
#'     minimum cv error (\code{p.1se}), and a data frame containing
#'     the loss for each fold, mean loss, and standard error of the
#'     loss for each value of p (\code{loss.df}).
#' 
#' @examples
#' data(treeda_example)
#' out.treedacv = treedacv(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     pvec = 1:10)
#' out.treedacv
#' @export
treedacv <- function(response, predictors, tree, folds = 5, pvec = 1:tree$Nnode,
                     k = nclasses - 1, center = TRUE, scale = TRUE, class.names = NULL, ...) {
    out = list()
    class(out) = "treedacv"
    if(is.null(class.names)) {
        class.names = unique(response)
    }
    nclasses = length(class.names)
    ## first create a vector giving the partitioning of the data
    if(length(folds) == 1) {
        ## divide the data as evenly as possible into the correct
        ## number of partitions
        partition = rep(1:folds, length.out = length(response))
        partition = sample(partition, size = length(partition), replace = FALSE)
    } else if(length(folds) != length(response)) {
        stop("'folds' must be either a single number or be of the same length as 'response'")
    } else {
        partition = folds
    }
    A = makeDescendantMatrix(tree)
    ## perform the cross validation
    loss.matrix = matrix(NA, nrow = length(pvec), ncol = max(partition) + 1)
    colnames(loss.matrix) = c(paste("Fold", 1:max(partition), sep = "."), "p")
    for(i in 1:max(partition)) {
        for(p.idx in 1:length(pvec)) {
            test.idx = which(partition == i)
            train.idx = which(partition != i)
            out.treeda = treeda(response[train.idx], predictors[train.idx,], tree,
                pvec[p.idx], k = k, center = center, scale = scale,
                check.consist = FALSE, class.names = class.names, A = A, ...)
            preds.treeda = predict.treeda(out.treeda, newdata = predictors[test.idx,],
                newresponse = response[test.idx], check.consist = FALSE)
            loss.matrix[p.idx,i] = mean(preds.treeda$classes != response[test.idx])
            loss.matrix[p.idx, ncol(loss.matrix)] = pvec[p.idx]
        }
    }

    ## summarize the results of cross validation
    loss.df = data.frame(loss.matrix)
    cvmeans = rowMeans(loss.matrix[,1:(ncol(loss.matrix) - 1)])
    cvses = apply(loss.matrix[,1:(ncol(loss.matrix) - 1)], 1,
        function(x) sd(x) / sqrt(length(x) - 1))
    min.idx = which.min(cvmeans)
    loss.df$means = cvmeans
    loss.df$ses = cvses
    out$p.min = pvec[min.idx]
    threshold = cvmeans[min.idx] + cvses[min.idx]
    out$p.1se = min(pvec[which(cvmeans <= threshold)])
    out$loss.df = loss.df
    return(out)   
}

#' Print treedacv objects
#' @param x \code{treedacv} object.
#' @param ... Not used
#' @method print treedacv
#' @export
print.treedacv <- function(x, ...) {
    cat("Output from cross-validation of treeda\n")
    cat("--------------------------------------\n")
    cat(paste("Value of p with minimum cv loss: ", x$p.min, "\n", sep = ""))
    cat(paste("Smallest p within 1 se of minimum cv loss: ", x$p.1se, "\n", sep = ""))
}


#' Plot a treedacv object
#'
#' Plots the cross-validation error with standard error bars.
#' @param x An object of class \code{treedacv}. 
#' @param ... Not used. 
#' @importFrom ggplot2 ggplot geom_point aes_string geom_errorbar ylab xlab
#' @examples
#' data(treeda_example)
#' out.treedacv = treedacv(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     pvec = 1:10)
#' plot(out.treedacv)
#' @export
#' @method plot treedacv
plot.treedacv <- function(x, ...) {
    p = ggplot(x$loss.df) + geom_point(aes_string(x = "p", y = "means")) +
        geom_errorbar(aes_string(x = "p", ymax = "means + ses", ymin = "means - ses"), width = .1) +
        ylab("Mean CV Error") + xlab("Sparsity")
    p
}
    

#' Make response matrix
#'
#' Create a dummy variable matrix for the response
#'
#' @param response A factor or character vector containing the
#' classes.
#' @param class.names A character vector giving the possible levels of
#' the factor. If NULL, it will be generated from the levels of
#' response.
#' @return A dummy variable matrix with column names giving the class
#' names.
#' @keywords internal
makeResponseMatrix <- function(response, class.names = NULL) {
    if(is.numeric(response)) {
        stop("'response' should be a factor or character vector.")
    }
    if(is.null(class.names)) {
        class.names = unique(response)
    }
    rm = matrix(NA, nrow = length(response), ncol = length(class.names))
    colnames(rm) = class.names
    for(j in 1:ncol(rm)) {
        rm[,j] = (response == colnames(rm)[j]) + 0
    }
    return(rm)
}


#' Make leaf coefficients
#'
#' For a set of coefficients defined on a matrix of (potentially
#' centered and scaled) leaf and node predictors, find the equivalent
#' set of coefficients on just the leaves.
#'
#' @param sda.out A fitted sda object
#' @param descendantMatrix A matrix describing the tree which was
#' used, element (i,j) indicates whether leaf i is a descendant of
#' node j.
#' @param means If the original predictor matrix was centered, the
#' means of the original predictor matrix, otherwise NULL.
#' @param sds If the original predictor matrix was scaled, the sds of
#' the original predictor matrix, otherwise NULL.
#' @return A list giving the coefficients on the leaves for each of
#' the discriminating axes and the intercepts for each of the
#' discriminating axes.
#' @importFrom Matrix Matrix
#' @keywords internal
makeLeafCoefficients <- function(sda.out, descendantMatrix, means, sds) {
    beta = sda.out$beta
    if(!is.null(sds)) {
        sdsSub = sds[sda.out$varIndex]
        beta = sweep(sda.out$beta, MARGIN = 1, STATS = sdsSub, FUN = "/")
    }
    if(is.null(means)) {
        intercepts = rep(0, ncol(beta))
    } else {
        intercepts = -means[sda.out$varIndex] %*% beta
    }
    leafCoef = Matrix(data = 0,nrow = ncol(descendantMatrix), ncol = ncol(beta))
    leafCoef[sda.out$varIndex,] = beta
    leafCoef = descendantMatrix %*% leafCoef
    return(list(beta = leafCoef, intercepts = intercepts))
}

#' Node coefficients to leaf coefficients
#'
#' General-purpose function for going from a coefficient vector on the
#' nodes to a coefficient vector on the leaves.
#'
#' @param coef.vec A vector containing coefficients on internal nodes plus leaves.
#' @param tree The phylogenetic tree.
#' @return A vector containing coefficients on the leaves. 
#' @export
nodeToLeafCoefficients <- function(coef.vec, tree) {
    descendantMatrix = makeDescendantMatrix(tree)
    leafCoef = descendantMatrix %*% coef.vec
    return(leafCoef)
}

#' Make descendant matrix
#'
#' Make a matrix describing the ancestry structure of a tree. Element
#' (i,j) indicates whether leaf i is a descendant of node j.
#'
#' @param tree A tree object of class phylo
#' @return A matrix describing the ancestry structure of a tree. 
#'
#' @importFrom Matrix Matrix
#' @keywords internal
makeDescendantMatrix <- function(tree) {
    etc = edgesToChildren(tree$edge)
    n = length(tree$tip.label)
    m = tree$Nnode
    A = Matrix::Matrix(data = 0, nrow = n, ncol = m + n)
    fillA = function(node, n) {
        if(node <= n) {
            A[node,node] <<- 1
        }
        else {
            children = etc[[node]]
            for(c in children) {
                fillA(c, n)
            }
            A[,node] <<- Matrix::rowSums(A[,children,drop=FALSE])
        }
    }
    fillA(n+1, n)
    return(A)
}

#' Make branch length vector
#'
#' Gets the branch lengths of the tree, with order the same as the
#' columns in makeDescendantMatrix.
#'
#' @param tree A tree object of class phylo
#' @return A vector of length ntips + nnodes, with ith element giving
#' the length of the branch above node i.
#'
#' @importFrom ape Ntip Nnode
#' @keywords internal
getBranchLengths <- function(tree) {
    branch_lengths = numeric(Ntip(tree) + Nnode(tree))
    branch_lengths[tree$edge[,2]] = tree$edge.length
    return(branch_lengths)
}


#' Make a matrix with predictors for each leaf and node
#'
#' Make a matrix with one predictor for each leaf and node in the
#' tree, where the node predictors are the sum of the leaf predictors
#' descending from them.
#'
#' @param leafPredictors A predictor matrix for the leaves: rows are
#'     samples, columns are leaves.
#' @param tree A phylogenetic tree describing the relationships
#'     between the species/leaves.
#' @return A predictor matrix for leaves and nodes together: rows are
#'     samples, columns are leaf/node predictors.
#' @export
makeNodeAndLeafPredictors <- function(leafPredictors, tree) {
    fullPredictors = as(leafPredictors, "matrix") %*% makeDescendantMatrix(tree)
    return(fullPredictors)
}

#' Print a treeda object
#' @param x \code{treeda} object.
#' @param ... Not used. 
#' @method print treeda
#' @export
print.treeda <- function(x, ...) {
    cat("An object of class treeda\n")
    cat("-------------------------\n")
    cat(paste(x$nPredictors,
              "predictors in the expanded space\nwere selected, corresponding to",
              x$nLeafPredictors, "\nleaves on the tree\n"))
    cat("-------------------------\n")
    cat("Confusion matrix:\n")
    with(x, print(table(truth = input$response, predicted = predictedClasses)))
}

#' Predict using new data
#'
#' Given a fitted \code{\link{treeda}} model, get the predicted
#' classes and projections onto the discriminating axes for new data.
#'
#' @param object Output from \code{\link{treeda}} function.
#' @param newdata New predictor matrix in the same format as the
#'     \code{predictor} argument to treeda. A matrix of predictor
#'     variables corresponding to the leaves of the tree and in the
#'     same order as the leaves of the tree.
#' @param newresponse New response vector, not required.
#' @param check.consist Check the consistency between the tree and
#'     predictor matrix?
#' @param ... Not used.
#' @return A list containing the projections of the new data onto the
#'     discriminating axes (\code{projections}), the predicted classes
#'     (\code{classes}), and the rss (\code{rss}, only included if the
#'     ground truth for the responses is available).
#' @importFrom mvtnorm dmvnorm
#' @examples
#' data(treeda_example)
#' out.treeda = treeda(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     p = 1)
#' ## Here we are predicting on the training data, in general this
#' ## would be done on a held out test set
#' preds = predict(out.treeda, newdata = treeda_example$predictors,
#'     newresponse = treeda_example$response)
#' ## make a confusion matrix
#' table(preds$classes, treeda_example$response)
#' @export
predict.treeda <- function(object, newdata, newresponse = NULL, check.consist = TRUE, ...) {
    if(check.consist) {
        checkPredictorsAndTree(newdata, object$input$tree)
    }
    out = list()
    out$projections = as(newdata %*% object$leafCoefficients$beta, "matrix") +
        matrix(1, nrow = nrow(newdata), ncol = 1) %*% object$leafCoefficients$intercept
    colnames(out$projections) = paste("Axis", 1:ncol(out$projections), sep = ".")
    ## get predicted classes
    nclasses = length(object$class.names)
    classLikelihood = matrix(0, nrow = nrow(out$projections),
        ncol = nrow(object$sda$theta))
    for(i in 1:nclasses) {
        cl = object$class.names[i]
        l = apply(out$projections, 1, function(x)
            with(object$classProperties,
                 dmvnorm(x, mean[[cl]], var[[cl]], log = TRUE)))
        classLikelihood[,i] = l + log(object$classProperties$prior[[cl]])
    }
    out$classes = object$class.names[apply(classLikelihood, 1, which.max)]
    ## compute the rss if ground truth for the responses is available
    if(!is.null(newresponse)) {
        responseMatrix = makeResponseMatrix(newresponse, object$class.names)
        out$rss = sum((responseMatrix %*% object$sda$theta - out$projections)^2)
    }
    return(out)
}


#' Coefficients from treeda fit
#'
#' Returns the coefficients from a treeda fit either in terms of the
#' leaves only or in terms of the nodes and leaves.
#' 
#' @param object An object of class \code{treeda}.
#' @param type Should the coefficients be in the leaf space or the
#'     node space?
#' @param ... Not used.
#'
#' @return A \code{\link[Matrix]{Matrix}} object containing the coefficients. 
#' 
#' @examples
#' data(treeda_example)
#' out.treeda = treeda(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     p = 1)
#' coef(out.treeda, type = "leaves")
#' coef(out.treeda, type = "nodes")
#' @importFrom Matrix Matrix
#' @export
coef.treeda <- function(object, type = c("leaves", "nodes"), ...) {
    type = match.arg(type)
    if(type == "leaves") {
        return(object$leafCoefficients$beta)
    }
    n = length(object$input$tree$tip.label)
    m = object$input$tree$Nnode
    beta = Matrix::Matrix(data = 0, nrow = n+m, ncol = ncol(object$leafCoefficients$beta))
    beta[object$sda$varIndex,] = object$sda$beta
    return(beta)
}

#' Compute properties of the classes
#'
#' For each class, computes the prior probabilities, means, and
#' variances for that class.
#'
#' @param response A vector containing the response for each observation.  
#' @param projections A matrix giving the projections of each
#' observation onto the discriminating axes.
#' @importFrom stats var
#' @keywords internal
makeClassProperties <- function(response, projections) {
    out = list()
    class.names = unique(response)
    out$prior = lapply(class.names, function(cl) mean(response == cl))
    out$mean = lapply(class.names, function(cl)
        colMeans(projections[response == cl,,drop = FALSE]))
    out$var = lapply(class.names, function(cl)
        var(projections[response == cl,,drop=FALSE]))
    names(out$prior) = names(out$mean) = names(out$var) = class.names
    return(out)
}
