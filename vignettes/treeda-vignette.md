---
title: "treeDA vignette"
author: "Julia Fukuyama"
date: "2017-07-24"
output:
  rmarkdown::html_document:
    toc: true
    toc_float: true
    theme: lumen
    keep_md: true
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{treeDA vignette}
  %\VignetteEncoding{UTF-8}
---




## Package overview

Here we will describe how to use the treeDA package. The package
provides functions to perform sparse discriminant analysis informed by
the tree. The method was developed for microbiome data, but it could
in principle be applied to any data with the same tree structure. The
idea behind the package is that when we have predictor variables which
are structured according to a tree, the mean values of the predictor
variables at each node in the tree are natural predictor variables,
and can be used in addition to the initial predictors defined at the
leaves. For microbiome data, this means using both the abundances of
the initial set of taxa as well as the abundances "pseudo-taxa", which
correspond to nodes in the tree and are the agglomeration of all the
taxa which descend from that node.

Without regularization, using both sets of predictors would yield an
ill-defined problem because the node predictors are linear
combinations of the leaf predictors. However, when we add
regularization, the problem becomes well-posed and we can obtain a
unique solution. Intuitively, the regularization allows us to
incorporate the intuition that a solution where one node is selected
is more parsimonious than one in which all the leaves descending from
that node are selected.

This package is based on the implementation of sparse discriminant
analysis implemented in the
[`sparseLDA`](https://cran.r-project.org/web/packages/sparseLDA/index.html)
package. The main function, `treeda`, creates the node and leaf
predictors, performs sparse discriminant analysis on the combination
of node and leaf predictors, and then translates the results back in
terms of leaf predictors only. The package also includes functions to
perform cross-validation and plotting, which will be demonstrated in
this vignette.

## Setup and data description

Our first step is to load the required packages and data. We will
illustrate the method on an antibiotic dataset (`AntibioticPhyloseq`)
provided by the package `adaptiveGPCA`. Note that no other elements of
the `adaptiveGPCA` package are used in this tutorial. The antibiotic
dataset consists of measurements taken from three subjects before,
during, and after taking each of two courses of an antibiotic. The
major groupings in the data are by subject (called `ind` in the
phyloseq object) and by the the antibiotic condition. The antibiotic
treatment is discretized into abx/no abx in a variable called `type`,
where abx corresponds to samples taken when the subject was taking the
antibiotic and the week following, and no abx corresponds to all the
other samples.

```r
library(treeDA)
library(ggplot2)
library(phyloseq)
library(adaptiveGPCA)
data(AntibioticPhyloseq)
theme_set(theme_bw())
```

## Model fitting

The main function in the package is called `treeda`. It takes a response
vector giving the classes to be separated, a matrix of predictor
variables which are related to each other by a tree, the tree which
describes the relationships between the predictor variables, and the
sparsity (p, the number of predictors to use). In the antibiotic
dataset, we have several potential discriminatory variables. One of
these describes whether the sample was taken during or immediately
after the subject was subjected to antibiotics, and we can try to find
taxa which discriminate between these two groups using the following
command:

```r
out.treeda = treeda(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 15)
```


## Model inspection
Here the output of the model is stored in an object called
`out.treeda`. The print function will give an overview of the fitted
model, including the number of predictors used and the confusion
matrix for the training data. 

```r
out.treeda
```

```
## An object of class treeda
## -------------------------
## 15 predictors in the expanded space
## were selected, corresponding to 903 
## leaves on the tree
## -------------------------
## Confusion matrix:
##         predicted
## truth    abx no abx
##   abx     58      9
##   no abx   7     88
```

From this, we see that 15 predictors were used (since this was what we
specified in the initial call to the function). These predictors
potentially include nodes in the tree (corresponding to taxonomic
clades) and leaves on the tree (corresponding to individual
species). The combination of nodes and leaves can be written purely in
terms of the leaves (or species, or OTUs), in which case the model is
using 903 of the leaves. This indicates that some of the nodes which
were selected as predictive were quite deep in the tree and
corresponded to large groups of taxa.

Finally, the confusion matrix shows us how well the model does on the
trainnig data: we see that a total of 16 cases were classified
incorrectly, split approximately evenly between cases which were
actually from the abx condition and those which were actually from the
no abx condition.


The object containing the output from the fit also contains other
information. These are:

- `means`: The mean value of each predictor. This is only included if
  the call to `treeda` included `center = TRUE`, in which case the
  means are stored so that new data can be centered using the mean
  values from the training data.

- `sds`: The standard deviation of each predictor. Like with the
  means, this is only included if the call to `treeda` included `scale
  = TRUE`, in which case the standard deviations are stored so that
  the new data can be scaled using the standard devaiations from the
  training data.

- `leafCoefficients`: A representation of the discriminating axis
using only the leaves. This is a list containing `beta`, which are the
coefficients, and `intercept`, which is the constant term.

- `input`: A list containing the response, predictors, and tree used
to fit the model.

- `nPredictors`: The number of predictors (in the node + leaf space)
used in the model.

- `nLeafPredictors`: The number of predictors in the leaf space used
in the model.

- `sda`: The sda object used in fitting the model.

- `class.names`: The names of the classes to be discriminated between.

- `projections`: The projections of the observations on the
discriminating axes.

- `classProperties`: The prior probabilities, mean in discriminating
space, and variance in the discriminating space of the classes.

- `predictedClasses`: Predicted classes for each observation.

- `rss`: Residual sum of squares: the sum of squared distances between
each observation and its class mean in the discriminating space.

## Sample plotting

Once we have fit the model, we can look at the samples projected onto
the discriminating axis. These projections are found in
`out.treeda$projections`, and we can see them plotted for the
antibiotic data below. In the figure below we also separate out the
samples by individual to see whether the model works better for some
individuals than others. We see that positive scores along the
discriminating axis correspond to the no abx condition, and that there
is some difference between the individuals but that the quality of the
model is approximately the same across the three subjects.


```r
ggplot(data.frame(sample_data(AntibioticPhyloseq), projections = out.treeda$projections)) +
    geom_point(aes(x = ind, y = projections, color = type))
```

![plot of chunk treeda-type-sample-plot](figure/treeda-type-sample-plot-1.png)

## Coefficient plotting
We can also look at the coefficient vector describing the
discriminating axis using the `plot_coefficients` function. This gives
a plot of the tree with the leaf coefficients aligned underneath. 

```r
plot_coefficients(out.treeda)
```

![plot of chunk treeda-type-coef-plot](figure/treeda-type-coef-plot-1.png)

For comparison, we can look at the results when we try to discriminate
between individuals instead of between the abx/no abx conditions. We
try this with the same amount of sparsity, p = 15. 

```r
out.treeda.ind = treeda(response = sample_data(AntibioticPhyloseq)$ind,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 15)
out.treeda.ind
```

```
## An object of class treeda
## -------------------------
## 30 predictors in the expanded space
## were selected, corresponding to 134 
## leaves on the tree
## -------------------------
## Confusion matrix:
##      predicted
## truth  D  E  F
##     D 56  0  0
##     E  0 52  0
##     F  0  0 54
```
In this case, since we have three classes we obtain two discriminating
axes, each of which uses 15 node or leaf predictors for a total of 30
predictors. This corresponds to only 85 leaves on the tree, indicating
that the nodes which were chosen corresponded to individual leaves or
to much smaller clades than when our aim was to discriminate between
the abx and no abx conditions. We can see this more clearly when we
look at the coefficient plot, where there are many more singleton
leaves with non-zero coefficients than we saw in the corresponding
plot for the abx/no abx model. Note that this model contains two
discriminating axes because we have three classes, while the abx/no
abx model had only one discriminating axis because there were two
classes. 

```r
plot_coefficients(out.treeda.ind)
```

![plot of chunk tereda-ind-coef-plot](figure/tereda-ind-coef-plot-1.png)

## Cross validation

We would often like to choose the sparsity level automatically instead
of manually. A common way of doing this is by cross validation, which
we have implemented in the function treedacv. It takes most of the
same arguments as as treeda: a vector containing the response, or the
classes for each of the observations, a matrix of predictors which are
related to each other by a tree, and the tree. In addition, the number
of folds for the cross validation needs to be specified (the folds
argument), and a vector giving the levels of sparsity to be compared
by cross validation (the pvec argument). The folds argument can be
given either as a single number, in which case the observations will
be partitioned into that number of folds, or as a vector assigning
each observation to a fold. In this case, the vector should have
length equal to the number of observations, and the elements in the
vector should be integers between 1 and the number of desired folds
assigning the observations to a fold.

Here we are using four-fold cross validation, discriminating between
the individuals in our dataset, and comparing levels of sparsity
between 1 and 15. When we print the output from `treedacv`, it tells us
both which value of p (amount of sparsity) corresponded to the minimum
CV error, and what the smallest value of p was which was within one
standard error of the minimum CV error. (The intuition behind using
this value of p instead of that with the minimum CV error is that we
would like the most parsimonious model which is statistically
indistinguishable from that with the minimum CV error). For us, the
minimum CV error is at 11, but if we were following the one standard
error rule we would use 7. 

```r
set.seed(0)
out.treedacv = treedacv(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq),
    folds = 4, pvec = 1:15)
out.treedacv
```

```
## Output from cross-validation of treeda
## --------------------------------------
## Value of p with minimum cv loss: 11
## Smallest p within 1 se of minimum cv loss: 7
```

The results from the cross validation are stored in
`out.treedacv$loss.df`. This data frame contains the CV error for each
fold, the mean CV error, and the standard error of the CV error for
each value of p. We can use this matrix to plot the CV error as a
function of the sparsity, or we can use the plotting function defined
by the package, as shown below.


```r
plot(out.treedacv)
```

![plot of chunk treeda-ind-plot-cv](figure/treeda-ind-plot-cv-1.png)

This plot confirms what we said earlier: 11 predictors corresponds to
the minimum cross validation error, and 7 predictors corresponds to
the sparsest solution which is within 1 standard error of the minimum
cross validation error.

We can then fit the model with 11 predictors to all the data and look
at the plot of the coefficients along the discriminating axis. 

```r
out.treeda.11 = treeda(response = sample_data(AntibioticPhyloseq)$type,
    predictors = otu_table(AntibioticPhyloseq),
    tree = phy_tree(AntibioticPhyloseq), p = 11)
out.treeda.11 
```

```
## An object of class treeda
## -------------------------
## 11 predictors in the expanded space
## were selected, corresponding to 892 
## leaves on the tree
## -------------------------
## Confusion matrix:
##         predicted
## truth    abx no abx
##   abx     58      9
##   no abx   8     87
```



```r
plot_coefficients(out.treeda.11)
```

![plot of chunk treeda-ind-plot-coef](figure/treeda-ind-plot-coef-1.png)

From the coefficient plot above, we might be interested in the
relatively large group of taxa with the largest positive
coefficients. Since the samples in the `abx` condition have positive
scores on the discriminating axis, taxa with positive coefficients are
over-represented in the `abx` condition. We can find out what these
are by examining the leaf coefficient vector. We first find all the
taxa which have the maximum positive coefficient, and then print out
the unique elements of the taxonomy table corresponding to those
taxa. We see that this is a group of 74 Lachnospiraceae. They are
mostly not annotated beyond the family level, but one is annotated as
being from the genus *Moryella*.


```r
taxa.max = which(out.treeda.11$leafCoefficients$beta == max(out.treeda.11$leafCoefficients$beta))
length(taxa.max)
```

```
## [1] 74
```

```r
unique(tax_table(AntibioticPhyloseq)[taxa.max,])
```

```
## Taxonomy Table:     [3 taxa by 8 taxonomic ranks]:
##          Domain     Phylum       Taxon_3        Taxon_4        
## Unc06b1n "Bacteria" "Firmicutes" "Clostridia_1" "Clostridiales"
## Unc072gp "Bacteria" "Firmicutes" "Clostridia_1" "Clostridiales"
## Unc05ycg "Bacteria" "Firmicutes" "Clostridia_1" "Clostridiales"
##          Taxon_5           Taxon_6          Taxon_7 Taxon_8
## Unc06b1n "Lachnospiraceae" "Incertae Sedis" ""      ""     
## Unc072gp "Lachnospiraceae" "uncultured"     ""      ""     
## Unc05ycg "Lachnospiraceae" "Moryella"       ""      ""
```
