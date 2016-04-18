#' Generate Correlation-Score Network Table from SEC protein.traces
#' @import data.table
#' @param rams protein.traces data.table with protein.traces, uniprot ids in column 1
#' @param cutoff Correlation score cutoff to apply before returning the output table
#' @param protcolnames Allows to set the column names of the protein-identifier containing columns.
#' 			Defaults to c("bait_id", "protein_id"), alternative, e.g.: protcolnames = c("uniprot_a", "uniprot_b").
#' @return data.frame with column headers as defined in protcolnames (uniprot_a uniprot_b) and
#'          cor_pearson and cor_spearman.
#' @export

generateCorNetwork <- function(protein.traces, cutoff = -1, protcolnames = c("uniprot_a", "uniprot_b")) {
  if (is.character(protein.traces)){
    protein.traces <- fread(protein.traces, header = TRUE)
  }

  data <- protein.traces[,2:ncol(protein.traces), with = FALSE]

  data <- as.matrix(data)
  rownames(data) <- protein.traces$protein_id

  corrmatrix_p <- cor(t(data))
  corrtable_p <- flattenCorrMatrix(corrmatrix_p)

  corrmatrix_s <- cor(t(data), method = "spearman")
  corrtable_s <- flattenCorrMatrix(corrmatrix_s)

  ## kendall deactivated, computation takes very long
  # corrmatrix_k <- cor(t(data), method = "kendall")
  # corrtable_k <- flattenCorrMatrix(corrmatrix_k)
  # corrtable <- cbind(corrtable_p, corrtable_s$cor, corrtable_k$cor)
  # names(corrtable)[5] <- "cor_kendall"

  corrtable <- cbind(corrtable_p, corrtable_s$cor)

  names(corrtable)[1] <- protcolnames[1]
  names(corrtable)[2] <- protcolnames[2]
  names(corrtable)[3] <- "cor_pearson"
  names(corrtable)[4] <- "cor_spearman"

  return(as.data.table(corrtable))
}
