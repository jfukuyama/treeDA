# treeDA Package

This package performs a version of sparse discriminant analysis which
takes into account the tree structure of the variables. It is intended
for use in microbiome data analysis, but in principle it could be used
in any situation in which the analyst has a number of classes to be
separated and a tree structure to the predictor variables.

The package is available on CRAN and can be installed using
```
install.packages("treeDA")
```
and then loaded with
```
library("treeDA")
```
Since `treeDA` was intended to be used to analyze microbiome data, it
depends on `phyloseq`, which is a bioconductor package and might need
to be installed separately. Installation of phyloseq is described on
the
[bioconductor website](http://bioconductor.org/packages/release/bioc/html/phyloseq.html),
and is accomplished with the following commands in R:
```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```

The main functions in the package are `treeda`, which fits a model
with a given sparsity level, and `treedacv`, which performs
cross-validation to choose the sparsity level. There is also a
function to plot the leaf coefficients obtained from the model aligned
to the tree (`plot_coefficients`). A full vignette showing an analysis
of microbiome data is available in the vignettes directory, or by
clicking
[here](http://htmlpreview.github.io/?https://github.com/jfukuyama/treeDA/blob/master/vignettes/treeda-vignette.html).

[![Travis-CI Build Status](https://travis-ci.org/jfukuyama/treeDA.svg?branch=master)](https://travis-ci.org/jfukuyama/treeDA)
