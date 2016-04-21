#' Benchmarks a Network against a Reference Network/Interaction set
#' @import data.table
#' @import ggplot2
#' @import ROCR
#' @param Network data.table or path to tab-separated .txt file or data.table with columns uniprot_a uniprot_b and confidence score
#' @param Reference Reference table, data.table or path to tab-separated text with headers containing info uniprot_a and uniprot_b for true interactions.
#'        Can also point to data.table with this information structure.
#' @param names_net Character vector specifiying the column names in Network containing the uniprot identifiers a and b and the confidence score.
#'        Defaults to c("uniprot_a", "uniprot_b")
#' @param scores_net Character vector of score column names in the network data. Defaults to the Bioplex scores c("pW", "pNI", "pInt").
#'
#' @param names_ref Character vector specifiying the column names in Network containing the uniprot identifiers a and b and the confidence score.
#'        Defaults to c("uniprot_a", "uniprot_b" , "score")
#' @param minscore_ref minimum score threshold in the reference.
#' @param n_shuffle_iter Number of shufflings to generate the negative test set.
#' @param output Character string, wether .pdfs and .csvs shall be written to working folder or whether the output should go to
#'        Rconsole. Options are "pdf_csv" or "Rconsole" with default "pdf_csv".
#' @return data.frame with column headers as defined by arg protcolnames and cor_pearson and cor_spearman
#' @export

benchmarkNetwork <- function(Network, Reference = "intact.txt",
                             names_net = c("UniprotA", "UniprotB"),
                             scores_net = c("pW", "pNI", "pInt"),
                             names_ref = c("uniprot_a", "uniprot_b", "miscore"),
                              minscore_ref = 0.5, n_shuffle_iter = 10,
                              output = "pdf_csv"){

  # input data treatment robust to supplying files or R objects
  if (is.character(Network)){
    net <- fread(Network)
    netname <- substr(Network, 1, nchar(Network)-4) #remove .txt
  } else {
    net <- as.data.table(Network)
    netname <- deparse(substitute(Network))
  }
  net <- subset(net, select = c(names_net, scores_net))
  setnames(net, names_net, c("uniprot_a", "uniprot_b"))

  if (is.character(Reference)){
    refset <- fread(Reference)
    refname <- substr(Reference, 1, nchar(Reference)-4) #remove .txt
    if (sum(grep("intact:", head(refset))) > 0){
      refset <- intactParser(refset)
    }
  } else {
    refset <- as.data.table(Reference)
    refname <- deparse(substitute(Reference))
  }
  ref <- subset(refset, select = names_ref)
  names(ref) <- c("uniprot_a", "uniprot_b", "score")
  ref$score <- as.numeric(ref$score)
  ref <- subset(ref, ref$score >= minscore_ref)

  # Create negative set
  # generate shuffled version of input data and label false
  net.s <- net
  for (i in 1:n_shuffle_iter){
    net.s <- transform(net.s, uniprot_a = sample(uniprot_a) )
    net.s <- transform(net.s, uniprot_b = sample(uniprot_b) )
  }
  negativeset <- unique(net.s)
  negativeset[, true_interaction:=FALSE]
  negativeset <- subset(negativeset, select = c("uniprot_a", "uniprot_b", "true_interaction"))

  # Create true set
  # Make Reference set bi-directional by appending reversed edges
  ref.bidirect <- rbind(ref[,c(1,2), with = FALSE], ref[,c(2,1), with = FALSE])
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
  for (i in 1:length(scores_net)){
    p <- ggplot(as.data.frame(trainingset), aes(color = true_interaction)) + aes_string(x = scores_net[i]) +
      geom_density() + theme_bw() + ggtitle(paste0(scores_net[i], "Score Density of forward and randomized Reference set"))
    plot(p)
  }
  if (output == "pdf_csv"){
    dev.off()
  }

  # Run ROCR/plot results
  # setkey(trainingset, cols = c("uniprot_a" , "uniprot_b"))
  trainingset <- as.data.frame(trainingset)
  if (output == "pdf_csv"){
    pdf(paste0("ROCplot_", netname, "_vs_", refname, ".pdf"))
  }

  auclist=c()
  fprlist001=c()
  cutofflist001=c()
  fprlist005=c()
  cutofflist005=c()

  i<-1
  while (i <= length(scores_net)) {
    pred <- prediction(trainingset[scores_net[i]],trainingset$true_interaction)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")

    if (i==1) {
      plot(perf, col=rainbow(length(scores_net))[i], main = paste0("ROC plot ", netname, " vs ", refname))
    }
    else {
      plot(perf, add = TRUE, col=rainbow(length(scores_net))[i])
    }

    auclist<-c(auclist,performance(pred, measure = "auc")@y.values[[1]])
    fprlist001<-c(fprlist001,perf@x.values[[1]][sum(perf@x.values[[1]]<=0.01)])
    cutofflist001<-c(cutofflist001,pred@cutoffs[[1]][sum(perf@x.values[[1]]<=0.01)])
    fprlist005<-c(fprlist005,perf@x.values[[1]][sum(perf@x.values[[1]]<=0.05)])
    cutofflist005<-c(cutofflist005,pred@cutoffs[[1]][sum(perf@x.values[[1]]<=0.05)])
    i<-i+1
  }
  abline( v = 0.01, lty = 3)
  abline( v = 0.05, lty = 4)
  legend("bottomright", scores_net, pch = 1, col=rainbow(length(scores_net)), cex=.75)
  auc<-data.frame("score"=scores_net,"AUC"=auclist,"FPR1"=fprlist001,"CutoffFPR1"=cutofflist001,"FPR5"=fprlist005,"CutoffFPR5"=cutofflist005)

  if (output == "pdf_csv"){
    dev.off()
  }

  table <- as.data.table((auc[with(auc, order(-AUC)), ]))
  # print(table(trainingset$true_interaction))
  if (output == "pdf_csv"){
    write.csv(table, file = paste0("ScoreCutoffs_", netname, "_vs_", refname, ".csv"), quote = FALSE, row.names = FALSE)
  }
  return(table)
}
