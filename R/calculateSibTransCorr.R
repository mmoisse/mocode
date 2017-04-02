#' Calculate sibling transition/fragment ion correlation in fragmention level
#' OpenSWATH/mProphet result table
#' @param data.transition Data.table of disaggregated OpenSWATH data An object of type \code{traces.obj}.
#' @param plot logical TRUE or FALSE 
#' @param PDF logical TRUE or FALSE 
#' @return data.transition Data.table of disaggregated OpenSWATH data including SibTransCorr column.
#' @export

calculateSibTransCorr <- function(data.transition,
                                plot = TRUE,
                                PDF = FALSE)
  {
  dataset.name <- deparse(substitute(data.transition))
  precursors = unique(data.transition$precursor_id)
  nprecursors <- length(precursors)
  SibTransCorr <- numeric(length = nprecursors)
  for (i in 1:nprecursors){
    
    message(paste("done", i, "of", nprecursors, "precursors"))
    indexpos <- precursors[i] == data.transition$precursor_id
    df <- data.transition[indexpos, .(FragmentIon, Intensity, Run)]
    df <- dcast(df, FragmentIon ~ Run, value.var = "Intensity")
    df.mat <- as.matrix(subset(df, select=-FragmentIon))
    rownames(df.mat) <- df$FragmentIon
    class(df.mat) <- 'numeric'
    df.mat_cor <- cor(t(df.mat))
    sibcorrs <- sapply(1:nrow(df.mat_cor), function(x){(sum(df.mat_cor[x,])-1)/(nrow(df.mat_cor)-1)})
    SibTransCorr[indexpos] <- sibcorrs
    
  }
  data.transition$SibTransCorr <- SibTransCorr
  
  # output plot
  if (plot){
    if (PDF){
      pdf(paste0(dataset.name, "_SibTransCorr_densityplot.pdf"))
    }
    plot(density(data.transition[grep("DECOY", ProteinName, invert = TRUE)]$SibTransCorr, na.rm = TRUE))
    lines(density(data.transition[grep("DECOY", ProteinName)]$SibTransCorr, na.rm = TRUE), lty = 2, col = "red")
    legend("topleft", legend = c("target fragment ions", "decoy fragment ions"), lty = c(1,2), col = c("black", "red"))
    if (PDF){
      dev.off()
    }
  }
  return(data.transition)
}



