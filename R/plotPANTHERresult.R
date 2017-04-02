# plot enrichment graphs
library(data.table)
library(ggplot2)
library(ggrepel)

plotPANTHERresult <- function(PANTHERresult = "panther_gocc.txt",
                              p_cutoff = 0.05){
  # define input
  if(class(PANTHERresult) == "character"){
    result <- fread(PANTHERresult)
  }
  
  # processing
  category <- names(result)[1]
  reflist <- names(result)[2]
  list <- names(result)[3]
  setnames(result, list, "in_list")
  setnames(result, category, "term")
  setnames(result, reflist, "in_reflist")
  result$category <- category 
  
  # clean up column names
  list <- strsplit(list, split = " ")[[1]][1]
  names(result)<- gsub(list, "", names(result))
  names(result) <- gsub(")", "", names(result))
  names(result) <- gsub(" \\(", "", names(result))
  names(result) <- gsub(" ", "_", names(result))
  names(result) <- gsub("\\-", "_", names(result))
  names(result) <- gsub("\\/", "_", names(result))
  
  # add list info columns
  result$list = list
  result$reflist = strsplit(reflist, split = " ")[[1]][1]
  
  # plotGoGraphs
  comparison_name <- paste0(gsub(".txt", "", unique(result$list)),
                            "_vs_",
                            gsub(".txt", "", unique(result$reflist)), 
                            "_", 
                            gsub(" ", "_", category))
  
  result <- result[term != "Unclassified (UNCLASSIFIED)"]
  result[, enrichment:=as.numeric(paste0(`over_under`, fold_Enrichment))]
  setorder(result, -enrichment)
  result <- result[!(is.na(enrichment))]
  result[, term_lean:=gsub(" \\(.*\\)", "", term)]$term_lean
  
  pdf(paste0(comparison_name, "_enrichment_scatter.pdf"), height = 5, width = 5)
  p = ggplot(result, aes(enrichment, -log10(P_value), label = term_lean, color = term_lean)) +
    geom_point(aes(size = in_list)) + geom_text_repel(data = result[P_value <= p_cutoff]) + theme_bw() + theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey") + 
    geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey")  +
    ggtitle(comparison_name)
  plot(p)
  dev.off()
  
  pdf(paste0(comparison_name, "_enrichment_bar.pdf"), height = 5, width = 5)
  q = ggplot(result, aes(x = reorder(term_lean, enrichment), enrichment, label = term_lean, fill = enrichment)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_continuous(low = "red", high = "blue") +
    geom_hline(yintercept = 0, linetype = 2) +
    xlab(unique(result$category)) +
    ylab("fold enrichment")
  plot(q)
  dev.off()
  
  return(list(table = result,
              scatterplot = p,
              barplot = q))
  
}


