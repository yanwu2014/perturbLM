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
#' @param pt.color Colors to use for data points
#' @return ggplot2 object with correlation plot
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
PlotCorrelation <- function(x, y, x.lab = NULL, y.lab = NULL, title = NULL, pts.use = NULL, use.label = F,
                            pts.label = NULL, font.size = 14, label.font.size = 4, pt.size = 1,
                            show.corr = F, box = T, alpha = 1, pt.color = NULL) {
  stopifnot(names(x) == names(y))

  if (is.null(pts.use)) {
    pts.use <- names(x)
  }
  corr.use <- paste("R = ", round(cor(x[pts.use], y[pts.use]), 2))

  if (!is.null(pt.color)) gg.df <- data.frame(x = x[pts.use], y = y[pts.use], color = pt.color[pts.use])
  else gg.df <- data.frame(x = x[pts.use], y = y[pts.use])
  gg.df$name <- names(x[pts.use])


  if (use.label) {
    gg.df$label <- pts.label
    gg.df$name[!pts.label] <- ""
  }

  if (show.corr) {
    if (is.null(title)) main.title <- corr.use
    else main.title <- paste(title, corr.use, sep = ": ")
  } else {
    main.title <- title
  }

  min.pt <- min(c(min(x[pts.use]), min(y[pts.use])))
  max.pt <- max(c(max(x[pts.use]), max(y[pts.use])))

  if (!is.null(pt.color)) ggobj <- ggplot(gg.df, aes(x, y, colour = color))
  else ggobj <- ggplot(gg.df, aes(x, y), color = "black")

  ggobj <- ggobj +
    geom_point(size = pt.size, alpha = alpha) +
    theme_classic() +
    theme(text = element_text(size = font.size), legend.title = element_blank()) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))

  if (box) {
    ggobj <- ggobj + xlim(c(min.pt, max.pt)) + ylim(c(min.pt, max.pt))
  }

  if (use.label) {
    ggobj <- ggobj + geom_text_repel(aes(x, y, label = name), size = label.font.size,
                                     colour = "black")
  }

  if (is.null(pt.color)) {
    ggobj <- ggobj + theme(legend.position = "none")
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



#' Violin plot (adapted from Seurat)
#'
#' @param feature Feature name
#' @param data Numeric vector of data to plot
#' @param cell.ident Factor with cell identities (must be same length as data)
#' @param do.sort Sort factor levels
#' @param y.max Set max expression level
#' @param size.x.use X axis text font size
#' @param size.y.use Y axis text font size
#' @param x.title X-axis title
#' @param size.title.use X & Y axis title font sizes
#' @param point.size.use Jitter point size
#' @param cols.use Color palette to use
#' @param y.log Log-scale gene expression
#' @param x.lab.rot Rotate x-axis labels
#' @param y.lab.rot Rotate y-axis labels
#' @param legend.position Legend position ("right" or "left")
#' @param remove.legend If true, hides legend
#'
#' @return A violin plot
#' @import ggplot2
#' @export
#'
PlotViolin <- function(feature, data, cell.ident, do.sort = F, y.max = NULL, size.x.use = 10,
                       size.y.use = 10, x.title = "Group", size.title.use = 11, point.size.use = 1,
                       cols.use = NULL, y.log = F, x.lab.rot = F, y.lab.rot = F, legend.position = "right",
                       remove.legend = F) {
  feature.name <- feature

  data <- data.frame(data)
  colnames(data) <- "feature"
  feature <- "feature"
  set.seed(seed = 42)
  data$ident <- cell.ident
  if (do.sort) {
    data$ident <- factor(data$ident, levels = names(rev(sort(tapply(data[, feature], data$ident, mean)))))
  }
  if (y.log) {
    noise <- rnorm(n = length(data[, feature]))/200
    data[, feature] <- data[, feature] + 1
  }
  else {
    noise <- rnorm(n = length(data[, feature]))/1e+05
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  }
  else {
    data[, feature] <- data[, feature] + noise
  }

  if(is.null(y.max)) {
    y.max <- max(data[, feature])
  }

  plot <- ggplot(data = data, mapping = aes(x = factor(ident), y = feature)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE, mapping = aes(fill = factor(ident))) +
    guides(fill = guide_legend(title = NULL)) + geom_jitter(height = 0, size = point.size.use) +
    xlab(x.title) + theme_classic() + ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"), legend.position = legend.position,
          axis.title.x = element_text(face = "bold", colour = "#990000", size = size.y.use))

  plot <- plot + ggtitle(feature.name)
  if (y.log) {
    plot <- plot + scale_y_log10()
  }
  else {
    plot <- plot + ylim(min(data[, feature]), y.max)
  }
  if (y.log) {
    plot <- plot + ylab(label = "Log Expression level")
  }
  else {
    plot <- plot + ylab(label = "Expression level")
  }
  if (!is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = size.x.use))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, size = size.y.use))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}



#' Dot plot visualization
#'
#' Intuitive way of visualizing geneset enrichment across different
#' groups (genotypes, clusters). The size of the dot encodes the percentage of
#' gsea enrichment within a group, while the color encodes the Fisher enrichment
#' cells within a group (blue is high).
#'
#' @param data.plot.df Dataframe with Geneset, Group, and enrichment info
#' @param color.col data.plot.df column to determine dot color
#' @param size.col data.plot.df column to determine dot size
#' @param row.col Plot rows (usually genesets)
#' @param col.col Plot columns (usually groups)
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
#' @param font.size x and y axis font sizes
#'
#' @return returns a ggplot2 object
#'
#' @import ggplot2
#' @export
#'
GenesetDotPlot <- function(data.plot.df, color.col, size.col, row.col = "Geneset", col.col = "Group",
                           cols.use = c("lightgrey", "blue"), col.max = 4, dot.scale = 4, scale.by = 'radius',
                           scale.min = NA, scale.max = NA, plot.legend = F, font.size = 12) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )

  data.plot <- data.frame(Color = data.plot.df[[color.col]], Size = data.plot.df[[size.col]],
                          Row = data.plot.df[[row.col]], Column = data.plot.df[[col.col]])
  data.plot$Color[data.plot$Color > col.max] <- col.max

  p <- ggplot(data = data.plot, mapping = aes(x = Column, y = Row)) +
    geom_point(mapping = aes(size = Size, color = Color)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = font.size, angle = 90),
          axis.text.y = element_text(size = font.size))

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



#' Make simple barplot from vector
#'
#' @param x Named numeric vector
#' @param y.lim Max y-axis
#' @param fill.color Bar color
#'
#' @return barplot
#' @import ggplot2
#' @export
#'
ggBarplot <- function(x, y.lim = NULL, fill.color = "lightgrey") {
  barplot.df <- data.frame(Y = x, X = factor(names(x), levels = names(x)),
                           color = fill.color)

  if(length(fill.color) == 1) {
    ggobj <- ggplot(data = barplot.df, aes(x = X, y = Y)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black", fill = fill.color)
  } else {
    ggobj <- ggplot(data = barplot.df, aes(x = X, y = Y, fill = color)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black")
  }


  ggobj <- ggobj + theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1, size = 14, angle = 90, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"))

  if(!is.null(y.lim)) {
    ggobj <- ggobj + coord_cartesian(ylim = y.lim)
  }
  ggobj
}
