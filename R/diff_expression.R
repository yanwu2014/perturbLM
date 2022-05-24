## Functions for computing single cell and pseudobulk differential expression

#' Compute per-cluster differential expression between conditions
#'
#' @param obj Seurat object
#' @param cluster.var Name of cluster metadata column
#' @param condition.var Name of condition metadata column
#' @param condition.1 Name of first condition to compare
#' @param condition.2 Name of second condition to compare
#' @param min.cells Minimum cells per cluster/condition combination
#' @param ... Additional keyword parameters to pass to FindMarkers
#'
#' @return Dataframe with per-cluster DE results
#' @export
ClusterWiseDE <- function(obj,
                          cluster.var,
                          condition.var,
                          condition.1,
                          condition.2,
                          min.cells = 10,
                          ...)
{
  require(Seurat)

  obj@meta.data[[cluster.var]] <- droplevels(as.factor(obj@meta.data[[cluster.var]]))
  obj@meta.data[[condition.var]] <- droplevels(as.factor(obj@meta.data[[condition.var]]))
  obj@meta.data$cluster_condition <- interaction(as.factor(obj@meta.data[[cluster.var]]),
                                                 as.factor(obj@meta.data[[condition.var]]),
                                                 sep = "||")
  meta.df <- obj@meta.data

  res_dict <- lapply(levels(meta.df[[cluster.var]]), function(cl) {
    cl.cond1 <- paste0(cl, "||", condition.1)
    cl.cond2 <- paste0(cl, "||", condition.2)

    cells.1 <- rownames(subset(meta.df, cluster_condition == cl.cond1))
    cells.2 <- rownames(subset(meta.df, cluster_condition == cl.cond2))

    if (length(cells.1) > min.cells && length(cells.2) > min.cells) {
      df <- FindMarkers(obj,
                        group.by = "cluster_condition",
                        ident.1 = cl.cond1,
                        ident.2 = cl.cond2,
                        ...)
      df$feature <- rownames(df)
      df$cluster <- cl
    } else {
      df <- NULL
    }

    return(df)
  })
  res_dict <- res_dict[!sapply(res_dict, is.null)]
  res_df <- do.call(rbind, res_dict)
  rownames(res_df) <- NULL

  return(res_df)
}
