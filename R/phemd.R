#' Visualize correlation between two vectors.
#' Adapted from PhEMD (https://bioconductor.org/packages/release/bioc/html/phemd.html)
#'
#' @param pheno.dict Named list of cell vectors for each genotype
#' @param clusters Cell type assignments
#' @param clusters.dists Matrix of distances between cell type centroids
#' @return Earth mover's distance between all genotypes
#'
#' @import pracma
#' @import transport
#' @export

computeEMD <- function(pheno.dict, clusters, clusters.dists) {
  require(swne)

  clusters <- factor(clusters)
  pheno.cluster.counts <- GenotypeClusterCounts(pheno.dict, UnflattenGroups(clusters))
  cluster_weights <- pheno.cluster.counts/rowSums(pheno.cluster.counts)
  clusters.dists <- clusters.dists[colnames(cluster_weights), colnames(cluster_weights)]

  # generate inhibitor distance matrix
  Y <- rep(0, (nrow(cluster_weights)-1)*nrow(cluster_weights)/2)
  counter <- 1
  for(i in seq_len(nrow(cluster_weights))){
    cur_f1_weights <- cluster_weights[i,]
    for(j in (i+1):nrow(cluster_weights)) {
      if(j > nrow(cluster_weights)) break #this doesn't automatically happen in R
      cur_f2_weights <- cluster_weights[j,]
      # Compute EMD
      flow <- transport(cur_f1_weights, cur_f2_weights, clusters.dists,
                        method='primaldual')
      curdist <- 0
      for(k in seq_len(nrow(flow))) {
        cur_penalty <- clusters.dists[flow[k,1], flow[k,2]]
        curdist <- curdist+cur_penalty*flow[k,3]
      }
      Y[counter] <- curdist
      counter <- counter + 1
    }
  }
  Y_sq <- squareform(Y)
  rownames(Y_sq) <- colnames(Y_sq) <- rownames(cluster_weights)
  return(Y_sq)
}
