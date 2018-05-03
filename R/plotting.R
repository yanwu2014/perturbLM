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
