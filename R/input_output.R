## Functions for reading/writing to file

#' Read in sample groups from a csv file format
#'
#' @param groups.file Name of groups file
#' @param sep.char Delimiter
#'
#' @return List of groups
#'
#' @export
#'
ReadGroups <- function(groups.file, sep.char = ",") {
  group.data <- read.table(groups.file, sep = sep.char, header = F, stringsAsFactors = F)
  groups <- group.data[[1]]
  full.groups.list <- lapply(rownames(group.data), function(i) sapply(strsplit(group.data[i,2], split = ',')[[1]], trimws))
  # full.groups.list <- lapply(full.groups.list, function(x) make.names(x))
  names(full.groups.list) <- groups
  return(full.groups.list)
}


#' Write sample groups from list to csv format
#'
#' @param groups.list List of groups
#' @param out.file Name of output file
#'
#' @export
#'
WriteGroups <- function(groups.list, out.file) {
  group.data <- sapply(groups.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  group.data <- sapply(names(group.data), function(x) paste(x, group.data[[x]], sep = ","))
  fileConn = file(out.file)
  writeLines(group.data, fileConn)
  close(fileConn)
}


#' Write sparse matrix using the 10X file format
#'
#' @param m Sparse matrix (features x barcodes)
#' @param out.dir Output directory
#' @param feature_id Vector of feature IDs (i.e. ENSEMBL IDs)
#' @param feature_short_name Vector of feature short names (i.e. gene symbols)
#'
#' @export
#'
Write10X <- function(m,
                     out.dir,
                     feature_id = NULL,
                     feature_short_name = NULL,
                     feature_type = "Gene Expression") {
  ## Create output directory
  dir.create(out.dir)

  ## Write barcodes
  write(colnames(m),
        file = paste0(out.dir, "/barcodes.tsv"),
        ncolumns = 1,
        sep = "\n")

  ## Write features
  feature_df <- data.frame(id = rownames(m),
                           short_name = rownames(m),
                           type = feature_type)
  if (!is.null(feature_id)) feature_df$id <- feature_id
  if (!is.null(feature_short_name)) feature_df$short_name <- feature_short_name
  write.table(feature_df,
              file = paste0(out.dir, "/features.tsv"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)

  ## Write matrix
  writeMM(m, file = paste0(out.dir, "/matrix.mtx"))

  ## Compress everything
  system(paste0("gzip ", out.dir, "/barcodes.tsv"))
  system(paste0("gzip ", out.dir, "/features.tsv"))
  system(paste0("gzip ", out.dir, "/matrix.mtx"))
}
