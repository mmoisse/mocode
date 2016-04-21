#' Benchmark Correlation-Score Network against a Reference interaction table
#' @import data.table
#' @import gridExtra
#' @import ggplot2
#' @import ROCR
#' @param corNetwork Path to tab-separated .txt file or data.table with column headers uniprot_a uniprot_b cor_pearson and cor_spearman
#' @param reference Reference table, tab-separated text with headers uniprot_a and uniprot_b for true interactions.
#'        Can also point to data.table with this information structure.
#' @param output Character string, wether .pdfs and .csvs shall be written to working folder or whether the output should go to
#'        Rconsole. Options are "pdf_csv" or "Rconsole" with default "pdf_csv".
#' @return data.frame with column headers as defined by arg protcolnames and cor_pearson and cor_spearman
#' @export

benchmarkCorNetwork <- function(corNetwork, reference = "corum.txt", output = "pdf_csv"){

  # input data treatment robust to supplying files or R objects
  if (is.character(corNetwork)){
    net <- fread(corNetwork)
    netname <- substr(corNetwork, 1, nchar(corNetwork)-4) #remove .txt
  } else {
    net <- as.data.table(corNetwork)
    netname <- deparse(substitute(corNetwork))
  }

  if (is.character(reference)){
    refset <- fread(reference)
    ref <- subset(refset, select = c("uniprot_a", "uniprot_b"))
    refname <- substr(reference, 1, nchar(reference)-4) #remove .txt
  } else {
    refset <- as.data.table(reference)
    refname <- deparse(substitute(reference))
    ref <- subset(refset, select = c("uniprot_a", "uniprot_b"))
  }


  # Create negative set
  # generate shuffled version of input data and label false
  net.s <- transform(net, uniprot_a = sample(uniprot_a) )
  net.s <- transform(net, uniprot_b = sample(uniprot_b) )
  negativeset <- unique(net.s)
  negativeset[, true_interaction:=FALSE]
  negativeset <- subset(negativeset, select = c("uniprot_a", "uniprot_b", "true_interaction"))

  # Create true set
  # Make reference set bi-directional by appending reversed edges
  ref.bidirect <- rbind(ref, ref[,c(2,1), with = FALSE])
  trueset <- unique(ref.bidirect)
  trueset[, true_interaction:=TRUE]
  trueset <- subset(trueset, select = c("uniprot_a", "uniprot_b", "true_interaction"))

  # merge into one trainingset
  truehits <- merge(net, trueset, by = c("uniprot_a", "uniprot_b"))
  negativehits <- merge(net, negativeset, by = c("uniprot_a", "uniprot_b"))
  trainingset <- rbind(truehits,negativehits[1:nrow(truehits)])

  # plot trainingset corr score density distributions
  if (output == "pdf_csv"){
    pdf(paste0("ReferenceSetScoresDensities_",netname, refname, ".pdf"), width = 12, height = 5.5)
  }
  pearson_dens <- ggplot(as.data.frame(trainingset), aes(x = cor_pearson, color = true_interaction)) +
  geom_density() + theme_bw()
  spearman_dens <- ggplot(as.data.frame(trainingset), aes(x = cor_spearman, color = true_interaction)) +
  geom_density()+ theme_bw()
  grid.arrange(pearson_dens, spearman_dens, ncol = 2, top = "Correlation coefficient density of training sets")
  if (output == "pdf_csv"){
    dev.off()
  }


  # Run ROCR
  pred.pearson <- prediction(trainingset$cor_pearson , trainingset$true_interaction)
  perf.pearson <- performance(pred.pearson, measure = "tpr", x.measure = "fpr")

  pred.spearman <- prediction(trainingset$cor_spearman , trainingset$true_interaction)
  perf.spearman <- performance(pred.spearman, measure = "tpr", x.measure = "fpr")

  # plot results
  if (output == "pdf_csv"){
    pdf(paste0("ROCplot_", netname, refname, ".pdf"))
  }
  plot(perf.pearson, col = "black", type = "l", main = paste("ROC curves of interaction recall of",netname, "vs", refname))
  plot(perf.spearman, add = TRUE, col = "red", type = "l")
  legend("bottomright", col = c("black", "red"), cex = 0.8,
         legend = c("cor_pearson", "cor_spearman"), lty = c(1,1))
  if (output == "pdf_csv"){
    dev.off()
  }

  #summarize and write out results underlying the curves
  result <- data.table(pearson_FPR = perf.pearson@x.values[[1]],
                      pearson_TPR = perf.pearson@y.values[[1]],
                      pearson_cutoff = perf.pearson@alpha.values[[1]],
                      spearman_FPR = perf.spearman@x.values[[1]],
                      spearman_TPR = perf.spearman@y.values[[1]],
                      spearman_cutoff = perf.spearman@alpha.values[[1]])
  message("some values may have been recycled in filling the result table")

    if (output == "pdf_csv"){
    write.csv(result, file = paste0("Performance_table_", netname, refname, ".csv"),
              quote = FALSE, row.names = FALSE)
  }
  return(result)
}
