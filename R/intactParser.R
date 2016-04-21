#' Intact Parser
#' @description Parses intact.txt or data.table to 3-column format with clean uniprot_a uniprot_b and miscore columns
#' @import data.table
#' @param intact intact data.table or path to .txt file
#' @param canonical whether uniprot entries shall be flattened to canonical, i.e. isoform flags (-x) removed
#' @return data.table with columns uniprot_a uniprot_b and miscore
#' @export
intactParser <- function(intact, canonical = TRUE){
  if (is.character(intact)){
    intact <- fread(intact)
  }
  intact.s <- subset(intact, select = c("#ID(s) interactor A", "ID(s) interactor B", "Confidence value(s)"))
  names(intact.s) <- c("uniprot_a", "uniprot_b", "miscore")
  intact.s <- unique(intact.s)
  intact.s <- intact.s[grep("uniprotkb:", intact.s$uniprot_a)]
  intact.s <- intact.s[grep("uniprotkb:", intact.s$uniprot_b)]
  intact.s[, uniprot_a:= gsub("uniprotkb:", "", uniprot_a)]
  intact.s[, uniprot_b:= gsub("uniprotkb:", "", uniprot_b)]
  intact.s[, miscore:= gsub("intact-miscore:", "", miscore)]

  if (canonical){
    intact.s[, uniprot_a:= sapply(uniprot_a, function(x){strsplit(x, split = "-")[[1]][1]})]
    intact.s[, uniprot_b:= sapply(uniprot_b, function(x){strsplit(x, split = "-")[[1]][1]})]
  }

  return(unique(intact.s))
}
