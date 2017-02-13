

#' Expand the background of a gtable.
#'
#' @param gtable A gtable object whose background needs to be expanded
#' to fill the whole space.
#'
#' @return A gtable object with a bigger background.

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
#' This method takes a ggplot of some data along the tips of the trees
#' and a ggplot of a tree and combines them. It assumes that you are
#' putting the tree on top and that the x axis for the plot has the
#' OTUs in the correct position (this can be found using the function
#' get_leaf_position).
#'
#' @param plot A plot of data about the OTUs with the x axis
#' corresponding to OTUs.
#' @param treePlot A plot of the tree.
#' @param treeHeight The relative amount of space in the plot the tree
#' should take up.
#' @return Prints the combined plot, returns a gtable object. 
#'
#' @importFrom ggplot2 ggplotGrob
#' @importFrom gtable gtable_add_rows
#' @importFrom gtable gtable_add_cols
#' @importFrom gtable gtable_add_grob
#' @importFrom grid unit
#' @export
combine_plot_and_tree = function(plot, tree.plot, tree.height = 5) {
    plot.grob = ggplotGrob(plot)
    tree.grob = ggplotGrob(tree.plot)
    tree.guide.idx = which(tree.grob$layout[,"name"] == "guide-box")
    tree.panel.idx = which(tree.grob$layout[,"name"] == "panel")
    plot.guide.idx = which(plot.grob$layout[,"name"] == "guide-box")
    plot.panel.idx = which(plot.grob$layout[,"name"] == "panel")
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
    gtable.bigger = gtable_add_rows(plot.grob, height = unit(tree.height, "null"), pos = 0)
    gtable.bigger.tree = gtable_add_grob(gtable.bigger, tree.panel.grob,
        t = 1, l = plot.panel.column, r = plot.panel.column)
    if(length(plot.guide.idx) == 0 & length(tree.guide.idx) == 1) {
        tree.guide.column = tree.guide.pos[,"l"]
        ## default position for gtable_add_columns is on the right,
        ## which is where we want it
        gtable.bigger.tree = gtable_add_cols(gtable.bigger.tree,
            width = tree.grob$widths[tree.guide.column])
        gtable.bigger.tree = gtable_add_grob(gtable.bigger.tree, tree.guide.grob,
            t = 1, l = dim(gtable.bigger.tree)[2], r = dim(gtable.bigger.tree)[2])
    } else if(length(plot.guide.idx) == 1 & length(tree.guide.idx) == 1) {
        gtable.bigger.tree = gtable_add_grob(gtable.bigger.tree, tree.guide.grob,
            t = 1, l = plot.guide.column, r = plot.guide.column)
    }
    ## add some padding around the top and right (right only if we
    ## needed to add an extra column for the tree guide)
    gtable.bigger.tree = gtable_add_rows(gtable.bigger.tree,
        height = plot.grob$heights[dim(plot.grob[1])], pos = 0)
    if(length(plot.guide.idx) == 0 & length(tree.guide.idx) == 1)
        gtable.bigger.tree = gtable_add_cols(gtable.bigger.tree, width = plot.grob$widths[1])
    gtable.bigger.tree = expand_background(gtable.bigger.tree)
    
    plot(gtable.bigger.tree)
    invisible(gtable.bigger.tree)
}


#' Get leaf positions from a tree layout
#'
#' Takes a tree, returns a vector with names describing the OTU and
#' entries giving the position of that OTU in the tree layout.
#'
#' @param tree A phylogenetic tree
#' @param ladderize FALSE for a non-ladderzied layout, TRUE or "right"
#' for a ladderized layout, "left" for a layout ladderized the other
#' way.
#' 
#' 
#' @export

get_leaf_position = function(tree, ladderize) {
    tree.layout = tree_layout(tree, ladderize = ladderize)
    dt = tree.layout$edgeDT
    dt.sub = subset(dt, !is.na(dt$OTU))
    out = data.frame(OTU = dt.sub$OTU, otu.pos = dt.sub$y)
    rownames(out) = out$OTU
    return(out)
}



#' Plots a tree and data
#'
#' Plots a tree on top and some data associated with the tips on the bottom.
#'
#' @param tree A tree of class phylo
#' @param data The data to plot, it needs to have row names that are
#' the same as the tip names of the tree.
#' @param annotation A vector of OTU names that should be marked.
#' @param ladderize How to lay out the tree, see tree_layout.
#' @param tree.height How big the tree is compared to the rest of the plots.
#' @param scales Passed to facet_grid, should be either "free_y" or
#' "fixed". 
#' @param annotation.color The color for the annotation dots.
#' @param barsep The distance of the taxonomy bar above the tree. 
#' @param barheight The width of the taxonomy bar. 
#' @param ... Additional arguments passed to plot_tree. 
#'
#' @return Plots the tree and data, returns a gtable object. 
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' 
#' @export
plot_tree_and_data <- function(tree, data, annotation = NULL, ladderize = TRUE,
                               tree.height = 5, scales = "fixed",
                               annotation.color = "firebrick1",
                               barsep = .02, barheight = .03, ...)
{
    
    combine_plot_and_tree(plot1, tree.plot, tree.height = tree.height)
}

