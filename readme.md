# treeDA Package

This package performs a version of sparse discriminant analysis which
takes into account the tree structure of the variables. It is intended
for use in microbiome data analysis, but in principle it could be used
in any situation in which the analyst has a number of classes to be
separated and a tree structure to the predictor variables.

The package can be installed using the command
```r
devtools::install_github("jfukuyam/treeDA")
```
and then loaded with
```r
library("treeDA")
```

The main functions in the package are `treeda`, which fits a model
with a given sparsity level, and `treedacv`, which performs
cross-validation to choose the sparsity level. There is also a
function to plot the leaf coefficients obtained from the model aligned
to the tree (`plot_coefficients`). A full vignette showing an analysis
of microbiome data is available in the vignettes directory.
