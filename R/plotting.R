#### Visualization functions ####


#' Visualize correlation between two vectors
#'
#' @param x X-axis vector
#' @param y Y-axis vector
#' @param x.lab X-axis title
#' @param y.lab Y-axis title
#' @param title Plot title
#' @param pts.use Data points to plot
#' @param labels Data points to label
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
#' @import ggrepel
#' @import ggpubr
#' @export
#'
PlotCorrelation <- function(x,
                            y,
                            x.lab = NULL,
                            y.lab = NULL,
                            title = NULL,
                            pts.use = NULL,
                            labels = NULL,
                            font.size = 14,
                            label.font.size = 4,
                            pt.size = 1,
                            x.lim = NULL,
                            y.lim = NULL,
                            show.corr = F,
                            alpha = 1,
                            pt.color = NULL)
{

  if (is.null(pts.use)) {
    pts.use <- intersect(names(x), names(y))
  }

  if (!is.null(pt.color)) {
    gg.df <- data.frame(x = x[pts.use], y = y[pts.use], color = pt.color[pts.use])
  }
  else {
    gg.df <- data.frame(x = x[pts.use], y = y[pts.use])
  }

  rownames(gg.df) <- pts.use
  gg.df$label <- ""
  gg.df[labels, "label"] <- labels

  if (!is.null(pt.color)) {
    ggobj <- ggplot(gg.df, aes(x, y, colour = color))
  }
  else {
    ggobj <- ggplot(gg.df, aes(x, y), color = "black")
  }

  ggobj <- ggobj +
    geom_point(size = pt.size, alpha = alpha) +
    theme_classic() +
    theme(text = element_text(size = font.size), legend.title = element_blank()) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))

  if (show.corr) {
    ggobj <- ggobj + stat_cor(method = "pearson")
  }

  if (!is.null(x.lim)) {
    ggobj <- ggobj + xlim(x.lim)
  }

  if (!is.null(y.lim)) {
    ggobj <- ggobj + xlim(y.lim)
  }

  if (!is.null(labels)) {
    ggobj <- ggobj + geom_text_repel(aes(x, y, label = label), size = label.font.size,
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





#' Boxplot of features
#'
#' @param obj Seurat object
#' @param feature Feature to plot
#' @param group.by Name of main grouping factor (the x-axis of boxplot)
#' @param fill.by Name of secondary grouping factor (fill color of the boxplot)
#' @param group.subset Plot only a subset of groups (default: NULL)
#' @param fill.subset Plot only a subset of fill groups (default: NULL)
#' @param pt.size Point size (set to zero to hide points)
#'
#' @return Boxplot of expression values
#' @import ggplot2
#' @export
PlotBox <- function(obj,
                    feature,
                    group.by,
                    fill.by,
                    group.subset = NULL,
                    fill.subset = NULL,
                    pt.size = 0.75)
{
  require(Seurat)
  df_plot <- obj@meta.data
  df_plot$expr <- GetAssayData(obj, "data")[feature,]
  df_plot[["group"]] <- df_plot[[group.by]]
  df_plot[["fill"]] <- df_plot[[fill.by]]

  if (!is.null(fill.subset)) {
    df_plot <- subset(df_plot, fill %in% fill.subset)
    df_plot[["fill"]] <- droplevels(df_plot[["fill"]])
  }

  if (!is.null(group.subset)) {
    df_plot <- subset(df_plot, group %in% fill.subset)
    df_plot[["group"]] <- droplevels(df_plot[["group"]])
  }

  ggplot(data = df_plot, mapping=aes(x=group, y=expr, fill=fill)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width=0.2),
               size = pt.size) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


#' Plot cluster compositions
#'
#' @param group.list Perturbation dictionary (in list format)
#' @param clusters Cluster annotations
#' @param samples Samples (i.e. mice, patients, biological replicates, default: NULL)
#' @param stabilize.var Apply variance stabilizing transformation (default: TRUE)
#' @param ctrl.group Control perturbation to compute logFC against (for paired samples only, default: NULL)
#' @param groups.plot Subset of perturbations to plot (default: NULL)
#' @param clusters.plot Subset of clusters to plot (default: NULL)
#' @param pt.size Point size (set to zero to hide points)
#'
#' @return Boxplot of cluster compositions
#' @import reshape2
#' @import ggplot2
#' @export
PlotClusterComp <- function(group.list,
                            clusters,
                            samples = NULL,
                            stabilize.var = TRUE,
                            ctrl.group = NULL,
                            groups.plot = NULL,
                            clusters.plot = NULL,
                            pt.size = 0.75)
{
  group_cluster_frac <- GroupClusterComp(group.list, clusters, samples,
                                         stabilize.var=stabilize.var,
                                         ctrl.group=ctrl.group)
  df_plot <- reshape2::melt(group_cluster_frac, value.name = 'proportion')

  if (!is.null(samples)) {
    df_plot[['group']] <- as.factor(sapply(as.character(df_plot[['Var1']]),
                                           ExtractField, field = 1, delim = '_'))
    df_plot[['sample']] <- as.factor(sapply(as.character(df_plot[['Var1']]),
                                            ExtractField, field = 2, delim = '_'))
  } else {
    df_plot[['group']] <- as.factor(df_plot[['Var1']])
  }

  df_plot[['cluster']] <- as.factor(df_plot[['Var2']])

  if (!is.null(groups.plot)) {
    df_plot <- subset(df_plot, group %in% groups.plot)
  }

  if (!is.null(clusters.plot)) {
    df_plot <- subset(df_plot, cluster %in% clusters.plot)
  }

  p <- ggplot(data = df_plot, mapping=aes(x=cluster, y=proportion, fill=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width=0.2),
               size = pt.size) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  if (!is.null(ctrl.group)) {
    p <- p + ylab('log2FC') + geom_hline(yintercept=0, color='red')
  }

  return(p)
}


#' Plot predicted values vs residuals to assess data linearity and heteroskedasticity
#' This requires densifying the expression matrices so it might be wise to only plot a subset of cells/genes.
#'
#' @param x Design matrix + covariates matrix
#' @param y Expression response
#' @param mfit Fitted glmnet model
#' @param nonzero.only Only plot residuals vs predicted for nonzero values
#'
#' @return KDE density plot of residuals vs predicted
#' @import glmnet
#' @export
#'
PlotResiduals <- function(x, y, mfit, nonzero.only = T) {
  y.pred <- predict(mfit, newx = x)
  y.pred <- as(pred.counts[,,1], "dgCMatrix")
  res <- y - y.pred

  df_temp <- data.frame(y_pred = as.vector(as.matrix(y.pred)),
                        res = as.vector(as.matrix(res)),
                        y = as.vector(y))
  if (nonzero.only) {
    df_temp <- subset(df_temp, y > 0)
  }

  graphics::smoothScatter(x = df_temp$y_pred,
                          y = df_temp$res,
                          xlab = "Predicted",
                          ylab = "Residuals",
                          nbin = 512,
                          nrpoints = 500,
                          pch = 20,
                          cex = 0.2,
                          col = "red")
}
