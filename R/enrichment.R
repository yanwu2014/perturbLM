#### Functions for handling genesets ####

#' Read genesets from gmt file to list format
#'
#' @param file.name Geneset file name
#' @return List of genesets
#'
#' @export
#'
LoadGenesets <- function(file.name) {
  geneset.to.use <- as.list(readLines(file.name))
  geneset.to.use <- lapply(geneset.to.use, function (v) strsplit(v, '\t')[[1]])
  geneset.names <- unlist(lapply(geneset.to.use, function(x) x[[1]]))
  geneset.to.use <- lapply(geneset.to.use, function(v) v[3:length(v)])
  names(geneset.to.use) <- geneset.names
  names(geneset.to.use) <- sapply(names(geneset.to.use), function(x) gsub(" ", "_", x))
  return(geneset.to.use)
}


#' Filter genesets by minimum, maximum number of genes
#'
#' @param genesets List of genesets
#' @param gene.name Gene names to filter
#' @param min.size Minimum genes per geneset
#' @param max.size Maxiumum genes per geneset
#' @return List of filtered genesets
#'
#' @export
#'
FilterGenesets <- function(genesets, gene.names, min.size = 5, max.size = 500) {
  genesets <- lapply(genesets, function (x) return(x[x %in% gene.names]))
  genesets <- .clean_genesets(genesets, min.size = min.size, max.size = max.size)
  return(genesets)
}


## Helper function for filtering genesets
.clean_genesets <- function(genesets, min.size = 5, max.size = 500, annot = FALSE) {
  genesets <- as.list(genesets)
  size <- unlist(lapply(genesets, length))
  genesets <- genesets[size > min.size & size < max.size]
  return(genesets)
}


#' Write genesets from list to gmt file
#'
#' @param genesets List of genesets
#' @param file.name GMT file name
#' @return Writes genesets to file.name
#'
#' @export
#'
WriteGenesets <- function(genesets, file.name) {
  genesets <- lapply(names(genesets), function(name) {x <- genesets[[name]]; x <- c(name,name,x); return(x);})
  n.cols <- 1.5*max(unlist(lapply(genesets, length)))
  empty <- lapply(genesets, write, file.name, append=TRUE, ncolumns = n.cols, sep = '\t')
}


#### Enrichment functions ####

#' Use Fisher's exact test to calculate enrichment for count data
#'
#' @param marker.genes Top marker genes
#' @param genesets List of genesets
#' @param n.background Number of background genes
#' @param correct FDR correction
#'
#' @return Dataframe with genesets, Fisher enrichment p-values (p.val), number of genes in the geneset (genes_in_set),
#'         size of each geneset (geneset_size)
#'
#' @export
#'
FisherEnrich <- function(marker.genes, genesets, n.background, correct = F) {
  results.df <- data.frame(genesets = names(genesets), p.val = numeric(length(genesets)),
                           genes_in_set = numeric(length(genesets)), geneset_size = numeric(length(genesets)),
                           stringsAsFactors = F)

  for (i in 1:length(genesets)) {
    set.genes <- genesets[[i]]
    set.size <- length(set.genes)
    in.set <- length(intersect(set.genes, marker.genes)) - 1
    out.set <- length(marker.genes) - in.set

    if (in.set < 0) { in.set <- 0 }
    cont.table <- matrix(c(in.set, set.size, out.set, n.background - set.size),
                         nrow = 2, byrow = T)

    results.df[i, "p.val"] <- fisher.test(cont.table, alternative = "greater")$p.value
    results.df[i, "geneset_size"] <- set.size
    results.df[i, "genes_in_set"] <- in.set
  }

  if (correct) {
    results.df$FDR <- p.adjust(results.df$p.val, method = "BH")
  }

  results.df <- results.df[order(results.df$p.val),]
  return(results.df)
}


#' Count over-enrichment for a list of gene markers
#'
#' @param markers.list List of top marker genes
#' @param genesets List of genesets
#' @param correct FDR correction
#'
#' @return Dataframe with group, genesets, Fisher enrichment p-values (p.val),
#'         number of genes in the geneset (genes_in_set),
#'         size of each geneset (geneset_size)
#'
#' @export
#'
MultipleFisherEnrich <- function(markers.list, genesets, correct = F) {
  all.genes <- unique(unlist(markers.list, F, F))
  genesets <- FilterGenesets(genesets, all.genes, min.size = 1, max.size = 10000)

  enrich.list <- lapply(names(markers.list), function(g) {
    res.df <- FisherEnrich(markers.list[[g]], genesets, length(all.genes), correct = correct)
    res.df$Group <- g; rownames(res.df) <- NULL;
    res.df
  })
  enrich.df <- do.call("rbind", enrich.list)
  # enrich.df$FDR <- p.adjust(enrich.df$p.val, method = "BY")
  enrich.df[order(enrich.df$p.val),]
}


#' GSEA over-enrichment test using liger
#'
#' @param scores.list List of scores to test enrichment for
#' @param genesets List of genesets
#' @param n.rand Number of randomizations
#' @param n.cores Number of cores
#' @param power GSEA exponent
#'
#' @return Dataframe with GSEA enrichment
#' @import liger
#' @export
#'
MultipleGSEAEnrich <- function(scores.list, genesets, n.rand = 1000, n.cores = 1, power = 1) {
  require(liger)

  enrich.list <- lapply(names(scores.list), function(g) {
    res <- bulk.gsea(scores.list[[g]], genesets, mc.cores = n.cores, n.rand = n.rand, power = power, rank = F,
                     skip.qval.estimation = T)
    res$Group <- g; res$genesets <- rownames(res); res$q.val <- NULL;
    res
  })
  enrich.df <- do.call("rbind", enrich.list); rownames(enrich.df) <- NULL;
  enrich.df$FDR <- p.adjust(enrich.df$p.val, method = "BY")
  enrich.df[order(enrich.df$p.val),]
}


#' Calculate chi-squared test of difference between each row of a counts matrix and the specified control row
#'
#' @param counts.df Counts matrix
#' @param ntc Control row
#'
#' @return Vector of p-values for each row
#' @export
#'
GenotypeClusterChisq <- function(counts.df, ntc) {
  apply(counts.df, 1, function(x) {
    chisq.test(rbind(x, counts.df[ntc,]))$p.value
  })
}


#' Use Fisher's exact test to calculate enrichment for each genotype in each cluster
#'
#' @param counts.df Genotype x cluster counts matrix
#' @return Genotype x cluster p-values matrix
#'
#' @export
#'
GenotypeClusterPvals <- function(counts.df) {
  cluster.counts <- colSums(counts.df)
  genotype.counts <- rowSums(counts.df)
  n.cells <- sum(cluster.counts)

  p.vals <- matrix(1, nrow = nrow(counts.df), ncol = ncol(counts.df),
                   dimnames = list(rownames(counts.df), colnames(counts.df)))

  for(i in 1:nrow(counts.df)) {
    for(j in 1:ncol(counts.df)) {
      x <- counts.df[i,j] - 1
      if (x < 0) { x <- 0; }
      cont.table <- matrix(c(x, cluster.counts[[j]], genotype.counts[[i]] - x, n.cells), nrow = 2, byrow = T)
      p.vals[i,j] <- fisher.test(cont.table, alternative = "greater")$p.value
    }
  }

  return(p.vals)
}
