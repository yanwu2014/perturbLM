#### Visualization functions ####


#' Visualize correlation between two vectors
#'
#' @param x X-axis vector
#' @param y Y-axis vector
#' @param x.lab X-axis title
#' @param y.lab Y-axis title
#' @param title Plot title
#' @param pts.use Data points to plot
#' @param use.label Label data points
#' @param pts.label Data points to label
#' @param font.size Axis and plot title font size
#' @param label.font.size Size to use for data point labels
#' @param pt.size Data point size
#' @param show.corr Display pearson correlation in title
#' @param box Force x and y vectors to be on the same scale
#' @param alpha Data point transparency
#' @return ggplot2 object with correlation plot
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
PlotCorrelation <- function(x, y, x.lab = NULL, y.lab = NULL, title = NULL, pts.use = NULL, use.label = F,
                            pts.label = NULL, font.size = 14, label.font.size = 4, pt.size = 1,
                            show.corr = F, box = T, alpha = 1) {
  stopifnot(names(x) == names(y))
  if (is.null(pts.use)) {
    pts.use <- names(x)
  }
  corr.use <- paste("R = ", round(cor(x[pts.use], y[pts.use]), 2))

  gg.df <- data.frame(x[pts.use],y[pts.use])
  gg.df$magnitude <- abs(x[pts.use]) + abs(y[pts.use])
  gg.df$magnitude <- as.numeric(gg.df$magnitude)
  gg.df$name <- names(x[pts.use])


  if (use.label) {
    gg.df$label <- pts.label
    gg.df$name[!pts.label] <- ""
    colnames(gg.df) <- c("x", "y", 'magnitude', 'name', 'label')
  } else {
    colnames(gg.df) <- c("x", "y", 'magnitude', 'name')
  }

  if (show.corr) {
    main.title <- paste(title, corr.use, sep = ": ")
  } else {
    main.title <- title
  }

  min.pt <- min(c(min(x[pts.use]), min(y[pts.use])))
  max.pt <- max(c(max(x[pts.use]), max(y[pts.use])))

  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(aes(colour = magnitude), size = pt.size, alpha = alpha) +
    scale_colour_gradient(low = 'grey', high = 'blue') + theme_classic() +
    theme(legend.position="none", text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title)
  if (box) {
    ggobj <- ggobj + xlim(c(min.pt, max.pt)) + ylim(c(min.pt, max.pt))
  }
  if (use.label) {
    ggobj <- ggobj + geom_text_repel(aes(x, y, label = name), size = label.font.size)
  }
  return(ggobj)
}

#' Visualize correlation between two vectors with density plot.
#' Useful for visualizing large datasets.
#'
#' @param v X-axis vector
#' @param w Y-axis vector
#' @param n.bins Number of hexagonal bins
#' @return ggplot2 object with hexbin plot
#'
#' @import ggplot2
#' @export
#'
PlotHexBin <- function(v, w, n.bins = 100) {
  ii <- intersect(names(v), names(w))
  r <- round(cor(v[ii], w[ii]), 3)
  df <- data.frame(x = v[ii], y = w[ii])
  ggobj <- ggplot(df) + geom_hex(aes(x, y), bins = n.bins, alpha = 1) +
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
    scale_fill_gradient(low = "deepskyblue3", high = "tomato3") + ggtitle(paste("R =", r))
  ggobj
}


#' Dot plot visualization
#'
#' Intuitive way of visualizing geneset enrichment across different
#' groups (genotypes, clusters). The size of the dot encodes the percentage of
#' gsea enrichment within a group, while the color encodes the Fisher enrichment
#' cells within a group (blue is high).
#'
#' @param gsea.df GSEA results
#' @param fisher.df Fisher enrichment results
#' @param genesets.plot Genesets to plot
#' @param groups.plot Groups to plot
#' @param cols.use Colors to plot, can pass a single character giving the name of
#' a palette from \code{RColorBrewer::brewer.pal.info}
#' @param col.max Maximum color value to set
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#' @param plot.legend plots the legends
#'
#' @return returns a ggplot2 object
#'
#' @import ggplot2
#' @export
#'
GenesetDotPlot <- function(gsea.df, fisher.df, genesets.plot, groups.plot, cols.use = c("lightgrey", "blue"),
                           col.max = 4, dot.scale = 4, scale.by = 'radius', scale.min = NA, scale.max = NA,
                           plot.legend = F) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )

  rownames(gsea.df) <- interaction(gsea.df$Group, gsea.df$genesets)
  rownames(fisher.df) <- interaction(fisher.df$Group, fisher.df$genesets)

  genesets.plot <- intersect(genesets.plot, unique(gsea.df$genesets))
  genesets.plot <- intersect(genesets.plot, unique(fisher.df$genesets))

  groups.plot <- intersect(groups.plot, unique(gsea.df$Group))
  groups.plot <- intersect(groups.plot, unique(fisher.df$Group))

  gsea.df <- subset(gsea.df, genesets %in% genesets.plot & Group %in% groups.plot)
  fisher.df <- subset(fisher.df, genesets %in% genesets.plot & Group %in% groups.plot)
  fisher.df <- fisher.df[rownames(gsea.df),]

  data.to.plot <- data.frame(group = gsea.df$Group, geneset = gsea.df$genesets, gsea.lp = -log(gsea.df$p.val),
                             fisher.lp = -log(fisher.df$p.val))
  data.to.plot <- data.to.plot[order(as.character(data.to.plot$group)),]
  data.to.plot$fisher.lp[data.to.plot$fisher.lp > col.max] <- col.max

  p <- ggplot(data = data.to.plot, mapping = aes(x = group, y = geneset)) +
    geom_point(mapping = aes(size = gsea.lp, color = fisher.lp)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  if (length(cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}
