#' plotPeptideProfiles
#' @description plot all peptide SEC chromatograms per protein
#' @param Traces list with $peptide.traces data tabl with protein_id, peptide_id and ordered quant values in the remaining columns
#' @param select_proteins character vector of protein_ids
#' @import data.table
#' @import reshape
#' @import ggplot2
#' @export

plotPeptideProfiles <- function(Traces, select_proteins = "all", output = "pdf"){
  objname <- deparse(substitute(Traces))
  data <- as.data.table(Traces$peptide.traces)
  labels <- data[,1:2, with = FALSE]
  quantdata <- as.matrix(data[,3:ncol(data), with = FALSE])
  rownames(quantdata) <- data$peptide_id
  if (select_proteins == "all"){
    proteins <- unique(data$protein_id)
  } else {
    proteins <- data$protein_id[data$protein_id %in% select_proteins]
  }
  nproteins <- length(proteins)
  if (output == "pdf"){
    pdf(paste0(objname, "_petideProfiles.pdf"), width = 12)
  }
  for (i in 1:nproteins){
    message(paste("processing plot", i, "of", nproteins))
    indexpos <- proteins[i] == data$protein_id
    df <- quantdata[indexpos,]
    #if it is only one peptide, manually construct the long table for plotting
    if (is.null(nrow(df))){
      peptide_id <- labels[indexpos]$peptide_id
      df_long <- data.table(peptide_id = peptide_id,
                            SEC_fraction = names(df),
                            Intensity = df)
    } else{
      df_long <- melt(df, varnames = c("peptide_id", "SEC_fraction"))
      setnames(df_long, "value", "Intensity")
    }
    p <- ggplot(df_long, aes(x = SEC_fraction, y = Intensity, group = peptide_id, colour = peptide_id)) +
      geom_line() + theme_minimal() + theme(legend.position="bottom") + ggtitle(paste("peptide SEC chromatograms for",proteins[i]))
    print(p)
  }
  if (output == "pdf"){
    dev.off()
  }
}
  