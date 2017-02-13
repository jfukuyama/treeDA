## tree definition in http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf

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

#' Tree-based SDA
#' 
#' Performs sparse discriminant analysis on tree-structured data by
#' augmenting the predictor matrix with additional predictors
#' corresponding to the nodes.
#'
#' @param response A factor or character vector giving the class to be
#' predicted.
#' @param predictors A matrix of predictor variables corresponding to
#' the leaves of the tree and in the same order as the leaves of the
#' tree.
#' @param tree A tree of class phylo.
#' @param p The number of predictors to use. 
#' @param k The number of components to use.
#' @param ... Additional arguments to be passed to sda
#' @importFrom sparseLDA sda
#' @importFrom Matrix colMeans
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
    predictions = predict(out, newdata = predictors, newresponse = response,
        check.consist = FALSE)
    out$predictedClasses = predictions$classes
    out$predictionError = predictions$predictionError
    out$rss = predictions$rss
    return(out)
}

#' treeda cross validation
#'
#' Performs cross-validation of a treeda fit.
#'
#' @param response The classes to be predicted. 
#' @param predictors A matrix of predictors corresponding to the tips
#' of the tree.
#' @param tree A tree object of class phylo.
#' @param folds Either a single number corresponding to the number of
#' folds of cross-validation to perform or a vector of integers
#' ranging from 1 to the number of folds desired giving the partition
#' of the dataset.
#' @param pvec The values of p to use.
#' @param k The number of discriminating axes to keep.
#' @param center Center the predictors?
#' @param scale Scale the predictors?
#'
#' @return A list with the value of p with minimum cv error, the
#' minimum value of p with in 1 e of the minimum cv error, a matrix
#' with the loss for each fold and each value of p, and vectors with
#' the mean and se of the loss for each value of p.
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
    
    ## perform the cross validation
    loss.matrix = matrix(NA, nrow = length(pvec), ncol = max(partition) + 1)
    colnames(loss.matrix) = c(paste("Fold", 1:max(partition), sep = "."), "p")
    for(i in 1:max(partition)) {
        for(p.idx in 1:length(pvec)) {
            test.idx = which(partition == i)
            train.idx = which(partition != i)
            out.treeda = treeda(response[train.idx], predictors[train.idx,], tree,
                pvec[p.idx], k = k, center = center, scale = scale,
                check.consist = FALSE, class.names = class.names, ...)
            preds.treeda = predict(out.treeda, newdata = predictors[test.idx,],
                newresponse = response[test.idx], check.consist = FALSE)
            loss.matrix[p.idx,i] = mean(preds.treeda$classes != response[test.idx])
            loss.matrix[p.idx, ncol(loss.matrix)] = pvec[p.idx]
        }
    }

    ## summarize the results of cross validation
    cvmeans = rowMeans(loss.matrix[,1:(ncol(loss.matrix) - 1)])
    cvses = apply(loss.matrix[,1:(ncol(loss.matrix) - 1)], 1,
        function(x) sd(x) / sqrt(length(x) - 1))
    min.idx = which.min(cvmeans)
    out$p.min = pvec[min.idx]
    threshold = cvmeans[min.idx] + cvses[min.idx]
    out$p.1se = min(pvec[which(cvmeans < threshold)])
    out$loss.matrix = loss.matrix
    out$cvmeans = cvmeans
    out$cvses = cvses
    return(out)   
}

#' Print treedacv object
#' 
#' @export
print.treedacv <- function(obj) {
    cat("Output from cross-validation of treeda\n")
    cat("--------------------------------------\n")
    cat(paste("Value of p with minimum cv loss: ", obj$p.min, "\n", sep = ""))
    cat(paste("Smallest p within 1 se of minimum cv loss: ", obj$p.1se, "\n", sep = ""))
}


#' Plot a treedacv object
#'
#' Gives a plot of the cross-validation error with standard error bars
#' 
#' @export
plot.treedacv <- function(obj) {
    df = data.frame(mean = obj$cvmeans, se = obj$cvse,
        p = obj$loss.matrix[,"p"])
    p = ggplot(df, aes(x = p, y = mean)) + geom_point() +
        geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = .1) +
        ylab("Mean cv loss")
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


#' Make descendant matrix
#'
#' Make a matrix describing the ancestry structure of a tree. Element
#' (i,j) indicates whether leaf i is a descendant of node j.
#'
#' @param tree A tree object of class phylo
#' @return A matrix describing the ancestry structure of a tree. 
#'
#' @importFrom Matrix Matrix
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
            A[,node] <<- Matrix::rowSums(A[,children])
        }
    }
    fillA(n+1, n)
    return(A)
}

#' Print a treeda object
#' 
#' @export
print.treeda <- function(obj) {
    cat("An object of class treeda\n")
    cat("-------------------------\n")
    cat(paste(obj$nPredictors,
              "predictors in the expanded space\nwere selected, corresponding to",
              obj$nLeafPredictors, "\nleaves on the tree\n"))
    cat("-------------------------\n")
    cat("Confusion matrix:\n")
    with(obj, print(table(truth = input$response, predicted = predictedClasses)))
    
}


#' Plot treeda object
#'
#' Provides a plot of the projections of the samples onto the
#' discriminating axes
#'
#' @param obj The treeda object to plot
#' @param type Plot the variables, the samples, or both? 
#' @param axes The axes to plot the samples along. 
#' @return A ggplot object
#' @export
plot.treeda <- function(obj, type = c("both", "variables", "samples"), axes = c(1,2)) {
    type = match.arg(type)
    if(type %in% c("both", "variables")) {
        
    }
    if(type %in% c("both", "samples")) {
        df = data.frame(obj$projections, class = obj$input$response)
        if(ncol(obj$projections) == 1) {
            psamples = ggplot(df) + geom_point(aes(y = class, x = Axis.1, color = class))
        } else {
            psamples = ggplot(df) +
                geom_point(aes_string(x = paste("Axis", axes[1], sep = "."),
                                      y = paste("Axis", axes[2], sep = "."),
                                      color = "class"))
        }
    }
    psamples
}


#' Predict using new data
#'
#' Given a fitted treeda model, get the predicted classes and
#' projections onto the discriminating axes for new data.
#'
#' @param treeda Output from treeda function
#' @param newdata New data
#' @return A list containing the projections of the new data onto the
#' discriminating axes, the predicted classes, and the rss (if the
#' ground truth for the responses is available).
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
predict.treeda <- function(treeda, newdata, newresponse = NULL, check.consist = TRUE) {
    if(check.consist) {
        checkPredictorsAndTree(newdata, treeda$input$tree)
    }
    out = list()
    out$projections = as(newdata %*% treeda$leafCoefficients$beta, "matrix") +
        matrix(1, nrow = nrow(newdata), ncol = 1) %*% treeda$leafCoefficients$intercept
    colnames(out$projections) = paste("Axis", 1:ncol(out$projections), sep = ".")
    ## get predicted classes
    nclasses = length(treeda$class.names)
    classLikelihood = matrix(0, nrow = nrow(out$projections),
        ncol = nrow(treeda$sda$theta))
    for(i in 1:nclasses) {
        cl = treeda$class.names[i]
        l = apply(out$projections, 1, function(x)
            with(treeda$classProperties,
                 dmvnorm(x, mean[[cl]], var[[cl]], log = TRUE)))
        classLikelihood[,i] = l + log(treeda$classProperties$prior[[cl]])
    }
    out$classes = treeda$class.names[apply(classLikelihood, 1, which.max)]
    ## compute the rss if ground truth for the responses is available
    if(!is.null(newresponse)) {
        responseMatrix = makeResponseMatrix(newresponse, treeda$class.names)
        out$rss = sum((responseMatrix %*% treeda$sda$theta - out$projections)^2)
    }
    return(out)
}


#' Get the coefficients from a treeda fit
#'
#' @export
coef.treeda <- function(obj, type = c("leaves", "nodes")) {
    type = match.arg(type)
    if(type == "leaves") {
        return(obj$leafCoefficients$beta)
    }
    n = length(obj$tree$tip.label)
    m = obj$tree$Nnode
    beta = Matrix::Matrix(data = 0, nrow = n+m, ncol = ncol(obj$leafCoefficients$beta))
    beta[obj$sda$varIndex,] = obj$sda$beta
    return(beta)
}


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
