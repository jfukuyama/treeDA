#' Expand the background of a gtable.
#'
#' @param gtable A gtable object whose background needs to be expanded
#' to fill the whole space.
#'
#' @return A gtable object with a bigger background.
#' @keywords internal
expand_background = function(gtable) {
    idx = which(gtable$layout[,"name"] == "background")
    if(length(idx) > 1) {
        warning("don't know what to do with a gtable with two backgrounds, returning the unmodified gtable")
        return(gtable)
    }
    if(length(idx) == 0) {
        warning("there is no background to expand, returning the unmodified gtable")
        return(gtable)
    }
    table.dim = dim(gtable)
    gtable$layout[idx, c("t", "l")] = c(1,1)
    gtable$layout[idx, "r"] = table.dim[2]
    gtable$layout[idx, "b"] = table.dim[1]
    return(gtable)
}

#' Method for combining two ggplots
#'
#' This method takes a ggplot of some data along the tips of the tree
#' and a ggplot of a tree and combines them. It assumes that you are
#' putting the tree on top and that the x axis for the plot has the
#' leaves in the correct position (this can be found using the
#' function \code{\link{get_leaf_position}}).
#'
#' @param plot A plot of data about the leaves with the x axis
#' corresponding to leaves.
#' @param tree.plot A plot of the tree.
#' @param tree.height The relative amount of space in the plot the tree
#' should take up.
#' @param print If true, the function will print the combined plot to
#' a graphics device, otherwise it will just return the gtable object
#' without printing.
#' @return Returns a \code{gtable} object. 
#'
#' @importFrom ggplot2 ggplotGrob
#' @importFrom gtable gtable_add_rows
#' @importFrom gtable gtable_add_cols
#' @importFrom gtable gtable_add_grob
#' @importFrom grid unit
#' @export
combine_plot_and_tree = function(plot, tree.plot, tree.height = 5, print = TRUE) {
    plot.grob = ggplotGrob(plot)
    tree.grob = ggplotGrob(tree.plot)
    tree.guide.idx = which(substr(tree.grob$layout[,"name"], 1, 9) == "guide-box")
    tree.panel.idx = which(substr(tree.grob$layout[,"name"], 1, 5) == "panel")
    plot.guide.idx = which(substr(plot.grob$layout[,"name"], 1, 9) == "guide-box")
    plot.panel.idx = which(substr(plot.grob$layout[,"name"], 1, 5) == "panel")
    if(length(tree.panel.idx) == 0) {
        stop("The tree plot isn't in the right form, don't know how to combine it.")
    }
    if(length(tree.guide.idx) == 1) {
        tree.guide.pos = tree.grob$layout[tree.guide.idx, 1:4, drop = FALSE]
        tree.guide.grob = tree.grob[tree.guide.pos[,"t"]:tree.guide.pos[,"b"], tree.guide.pos[,"l"]:tree.guide.pos[,"r"]]
    }
    if(length(tree.panel.idx) == 1) {
        tree.panel.pos = tree.grob$layout[tree.panel.idx, 1:4, drop = FALSE]
        tree.panel.grob = tree.grob[tree.panel.pos[,"t"]:tree.panel.pos[,"b"], tree.panel.pos[,"l"]:tree.panel.pos[,"r"]]
    }
    if(length(plot.guide.idx) > 0) {
        ## I'm assuming the guide takes only one cell, probably this
        ## is always true
        plot.guide.column = plot.grob$layout[plot.guide.idx, "l"]
    }
    if(length(plot.panel.idx) > 0) {
        plot.panel.column = unique(plot.grob$layout[plot.panel.idx, "l"])
        if(length(plot.panel.column) > 1)
            stop ("The data plot has too many data panels, don't know how to combine it with the tree plot.")
    }
    ## strip the white space on top here
    gtable.bigger = gtable_add_rows(plot.grob, heights = unit(tree.height, "null"), pos = 0)
    gtable.bigger.tree = gtable_add_grob(gtable.bigger, tree.panel.grob,
        t = 1, l = plot.panel.column, r = plot.panel.column)
    if(length(plot.guide.idx) == 0 & length(tree.guide.idx) == 1) {
        tree.guide.column = tree.guide.pos[,"l"]
        ## default position for gtable_add_columns is on the right,
        ## which is where we want it
        gtable.bigger.tree = gtable_add_cols(gtable.bigger.tree,
            widths = tree.grob$widths[tree.guide.column])
        gtable.bigger.tree = gtable_add_grob(gtable.bigger.tree, tree.guide.grob,
            t = 1, l = dim(gtable.bigger.tree)[2], r = dim(gtable.bigger.tree)[2])
    } else if(length(plot.guide.idx) == 1 & length(tree.guide.idx) == 1) {
        gtable.bigger.tree = gtable_add_grob(gtable.bigger.tree, tree.guide.grob,
            t = 1, l = plot.guide.column, r = plot.guide.column)
    }
    ## add some padding around the top and right (right only if we
    ## needed to add an extra column for the tree guide)
    gtable.bigger.tree = gtable_add_rows(gtable.bigger.tree,
        heights = plot.grob$heights[dim(plot.grob[1])], pos = 0)
    if(length(plot.guide.idx) == 0 & length(tree.guide.idx) == 1)
        gtable.bigger.tree = gtable_add_cols(gtable.bigger.tree, widths = plot.grob$widths[1])
    gtable.bigger.tree = expand_background(gtable.bigger.tree)
    if(print) {
        plot(gtable.bigger.tree)
        invisible(gtable.bigger.tree)
    }
    return(gtable.bigger.tree)
}


#' Get leaf positions from a tree layout
#'
#' Takes a tree, returns a vector with names describing the leaves and
#' entries giving the position of that leaf in the tree layout.
#'
#' @param tree A tree of class \code{phylo}. 
#' @param ladderize FALSE for a non-ladderzied layout, TRUE or "right"
#' for a ladderized layout, "left" for a layout ladderized the other
#' way.
#' @importFrom phyloseq tree_layout
#' @export
get_leaf_position = function(tree, ladderize) {
    tree.layout = tree_layout(tree, ladderize = ladderize)
    dt = tree.layout$edgeDT
    dt.sub = subset(dt, !is.na(dt$OTU))
    out = data.frame(OTU = dt.sub$OTU, otu.pos = dt.sub$y)
    rownames(out) = out$OTU
    return(out)
}


#' Plot the discriminating axes from treeda
#'
#' Plots the leaf coefficients for the discriminating axes in a fitted
#' \code{treeda} model aligned under the tree. 
#'
#' @param out.treeda The object resulting from a call to
#'     \code{\link{treeda}}.
#' @param remove.bl A logical, \code{TRUE} if the tree should be plotted
#'     after setting all branch lengths equal to the same value or
#'     not. The plots tend to look nicer when all the branch lengths
#'     are the same, and the branch length information is not used in
#'     the model.
#' @param ladderize Layout parameter for the tree.
#' @param tree.height The height of the tree relative to the height of
#'     the plot below.
#' @return A plot of the tree and the coefficients.
#' @importFrom grid grid.draw
#' @importFrom phyloseq plot_tree
#' @importFrom ggplot2 coord_flip scale_x_reverse facet_grid
#'     element_blank aes_string ylab theme
#' @examples
#' data(treeda_example)
#' out.treeda = treeda(response = treeda_example$response,
#'     predictors = treeda_example$predictors,
#'     tree = treeda_example$tree,
#'     p = 1)
#' plot_coefficients(out.treeda)
#' @export
plot_coefficients <- function(out.treeda, remove.bl = TRUE, ladderize = TRUE, tree.height = 2) {
    tr = out.treeda$input$tree
    if(remove.bl) {
        tr$edge.length = rep(1, length(tr$edge.length))
    }
    tree.plot = plot_tree(tr, ladderize = ladderize) + coord_flip() + scale_x_reverse()
    leaf.position = get_leaf_position(out.treeda$input$tree, ladderize = ladderize)$otu.pos
    coef = as(out.treeda$leafCoefficients$beta, "matrix")
    colnames(coef) = paste("Axis", 1:ncol(coef))
    df = data.frame(coef, leaf.position)
    df = reshape2::melt(df, id.vars = "leaf.position")
    coef.plot = ggplot(df) +
        geom_point(aes_string(x = "leaf.position", y = "value")) +
        facet_grid(variable ~ .) +
        ylab("Coefficient value") + 
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
    p = combine_plot_and_tree(coef.plot, tree.plot, tree.height = tree.height, print = FALSE)
    grid::grid.draw(p)
    invisible(p)
}
