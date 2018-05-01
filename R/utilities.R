#' Convert dataframe to matrix, specifying all column names
#'
#' @param df Dataframe
#' @param output.name Column of df to use as matrix values
#' @param row.col Column of df to use as matrix rows
#' @param col.col Column of df to use as matrix columns
#'
#' @return Matrix
#' @importFrom reshape2 acast
#' @export
#'
UnflattenDataframe <- function(df, output.name, row.col = 'Gene', col.col = 'Group') {
  df <- df[c(row.col, col.col, output.name)]
  colnames(df) <- c('row', 'column', output.name)
  mat.out <- reshape2::acast(df, row~column, value.var = output.name)
  return(mat.out)
}
UnflattenDataframe <- compiler::cmpfun(UnflattenDataframe)

