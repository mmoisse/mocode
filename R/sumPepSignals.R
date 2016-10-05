#' sumPepSignals
#' @param widepepdata wide peptide level quant data.table with label columns protein_id and peptide_id
#' @return Wide protein level table with peptide signals summed per protein_id
#' @export
sumPepSignals <- function(widepepdata){
  long <- melt(widepepdata)
  long.protein <- aggregate(long$value, by = list(long$protein_id, long$variable), FUN = sum)
  wide.protein <- dcast(long.protein, Group.1 ~ Group.2)
  names(wide.protein)[1]<-"protein_id"
  wide.protein
}
