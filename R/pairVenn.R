#' pairVenn
#'  @description Draws a pairwise Venn Diagram of two character vectors with overlap area
#'  corresponding to overlap size.
#'  @param a Character vector a
#'  @param b Character vector b
#'  @param output character selection either "R" or "pdf"
#'  @import VennDiagram
#'  @import gridExtra
#'  @export

pairVenn <- function(a, b, output = "R"){
  require(gridExtra)
  require(VennDiagram)
  namea <- deparse(substitute(a))
  nameb <- deparse(substitute(b))
  ua <- unique(a)
  ub <- unique(b)
  mainlabel = paste("Pairwise Venn diagram of list", namea, "and", nameb)
  
  pairvenn <- draw.pairwise.venn(area1 = length(ua),
                                 area2 = length(ub),
                                 cross.area = length(intersect(ua, ub)),
                                 category = c(namea, nameb),
                                 cat.col = c("blue", "red"),
                                 col = c("blue", "red"),
                                 fill = c("blue", "red"),
                                 alpha = 0.2,
                                 margin = rep(0.05, 4))
  
  if (output == "R"){
    grid.arrange(gTree(children=pairvenn), top=mainlabel)
  }
  
  if (output == "pdf"){
    pdf(paste(mainlabel, ".pdf"))
    grid.arrange(gTree(children=pairvenn), top=mainlabel)
    dev.off()
  }
}