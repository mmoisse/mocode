#' Flattens Correlation Matrix to table with columns row, column and cor
#' @param cormat correlation matrix with row- and colnames as produced by cor() function
#' @return data.frame with columns row, column and cor
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}
