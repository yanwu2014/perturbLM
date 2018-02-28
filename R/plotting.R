#### Visualization functions ####

#' Visualize t-SNE or PCA plots
#'
#' @param dim.x X-axis vector
#' @param dim.y Y-axis vector
#' @param clusters Cell clusters
#' @param x.lab X-axis title
#' @param y.lab Y-axis title
#' @param main.title Plot title
#' @param pt.size Data point size
#' @param font.size Font size for axis titles
#' @param alpha Data point transparency
#' @param do.label Whether or not to label the clusters on the plot
#' @param label.size Size of cluster labels
#' @param show.legend Display legend
#' @param label.points Label each individual points
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
PlotDims <- function(dim.x, dim.y, clusters = NULL, x.lab = "tsne1", y.lab = "tsne2", main.title = NULL,
                     pt.size = 1.0, font.size = 12, alpha = 1.0, do.label = T, label.size = 4,
                     show.legend = T, label.points = F) {
  gg.df <- data.frame(x = dim.x, y = dim.y, clusters = clusters, pts.label = names(dim.x))
  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(size = pt.size, alpha = alpha, aes(colour = clusters)) +
    theme_void() + theme(text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title)

  if (do.label && is.factor(clusters)) {
    group.pts.x <- tapply(gg.df$x, gg.df$clusters, function(v) {
      median(v)
    })
    group.pts.y <- tapply(gg.df$y, gg.df$clusters, function(v) {
      median(v)
    })
    group.pts <- data.frame(x = group.pts.x, y = group.pts.y)
    group.pts$ident <- levels(gg.df$clusters)

    ggobj <- ggobj + geom_point(data = group.pts, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
      geom_text_repel(data = group.pts, mapping = aes(label = ident), size = label.size)
  }

  if (!show.legend) {
    ggobj <- ggobj + theme(legend.position = "none")
  }

  if (label.points) {
    ggobj <- ggobj + geom_text_repel(mapping = aes(label = pts.label), size = label.size)
  }

  ggobj
}


#' Plot matrix as heatmap. Options for clustering rows/columns
#' @param m Input matrix to visualize
#' @param rescaling Rescale the matrix. Options include row/column wise scaling, or a custom scaling function
#' @param clustering Cluster rows/columns or both
#' @param labCol Label columns
#' @param labRow Label rows
#' @param border Plot heatmap border
#' @param heatscale Colormap to use
#' @param legend.title Colormap title
#' @param x.lab.size X-axis font size
#' @param y.labl.size Y-axis font size
#' @return ggplot2 object with heatmap plot
#'
#' @import ggplot2
#' @import reshape2
#' @export
#'
ggheat <- function(m, rescaling = 'none', clustering = 'none', labCol = T, labRow = T, border = F,
                   heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato'), legend.title = NULL,
                   x.lab.size = 8, y.lab.size = 8) {
  require(reshape)
  require(ggplot2)

  ## you can either scale by row or column not both!
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct

  if(is.function(rescaling)) {
    m = rescaling(m)
  } else
  {
    if(rescaling == 'column')
      m = scale(m, center=T)
    if(rescaling == 'row')
      m = t(scale(t(m),center=T))
  }

  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function

  if(is.function(clustering)) {
    m = clustering(m)
  } else {
    if(clustering=='row')
      m = m[hclust(dist(m))$order, ]
    if(clustering=='column')
      m = m[ ,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m = m[hclust(dist(m))$order, hclust(dist(t(m)))$order]
  }
  ## this is just reshaping into a ggplot format matrix and making a ggplot layer

  rows = dim(m)[1]
  cols = dim(m)[2]
  melt.m = cbind(rowInd=rep(1:rows, times = cols), colInd = rep(1:cols, each = rows), melt(m))
  g = ggplot(data = melt.m)

  ## add the heat tiles with or without a white border for clarity
  if(border == TRUE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value), colour = 'white')
  if(border == FALSE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value))

  ## add axis labels either supplied or from the colnames rownames of the matrix
  if(labCol == T)
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = colnames(m), expand = c(0.005,0))
  if(labCol == F)
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = rep('', cols))
  if(labRow == T)
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rownames(m), expand = c(0.005,0))
  if(labRow == F)
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rep('', rows))

  ## get rid of grey panel background and gridlines
  g2 = g2 + theme(panel.grid.minor = element_line(colour = NA), panel.grid.major = element_line(colour=NA),
                  panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_text(angle = 90, hjust = 1, size = x.lab.size),
                  axis.ticks = element_blank(), axis.text.y = element_text(size = y.lab.size))

  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return(g2 + scale_fill_gradient2(low = heatscale[1], mid = heatscale[2], high = heatscale[3], guide = guide_colorbar(title = legend.title)))
}


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
#' @param color.magnitude Color points by x + y magnitude
#' @return ggplot2 object with correlation plot
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#'
ggcorrelation <- function(x, y, x.lab = NULL, y.lab = NULL, title = NULL, pts.use = NULL, use.label = F,
                          pts.label = NULL, font.size = 14, label.font.size = 4, pt.size = 1,
                          show.corr = F, box = T, alpha = 1, color.magnitude = T) {
  stopifnot(names(x) == names(y))
  if (is.null(pts.use)) {
    pts.use <- names(x)
  }
  corr.use <- paste("R = ", round(cor(x[pts.use], y[pts.use]), 2))

  gg.df <- data.frame(x[pts.use],y[pts.use])
  gg.df$name <- names(x[pts.use])

  if (color.magnitude) {
    gg.df$magnitude <- abs(x[pts.use]) + abs(y[pts.use])
  } else {
    gg.df$magnitude <- 1
  }

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
    scale_colour_gradient(low = 'grey', high = 'darkblue') + theme_classic() +
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
gghexbin <- function(v, w, n.bins = 100) {
  ii <- intersect(names(v), names(w))
  r <- round(cor(v[ii], w[ii]), 3)
  df <- data.frame(x = v[ii], y = w[ii])
  ggobj <- ggplot(df) + geom_hex(aes(x, y), bins = n.bins, alpha = 1) +
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
    scale_fill_gradient(low = "skyblue", high = "tomato") + ggtitle(paste("R =", r))
  ggobj
}
