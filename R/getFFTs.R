#' Get fraction of false targets
#' @description Retrieve m/py-Prophet-model-estimated fraction of false targets (FFT) per run
#' @param folder folder where the *full_stat.csv files are
#' @param pattern regex pattern for file grabbing
#' @export
getFFTs <- function(folder = "pyprophet_stats/", pattern = "stats.csv"){
  statsfiles <- list.files(folder)[grep(pattern, list.files(folder))]
  qvalue1 <- numeric()
  for (i in 1:length(statsfiles)){
    file <- read.csv(paste0(folder, statsfiles[i]), sep = "\t")
    qvalue1[i] <- file$qvalue[1]
  }

  return(data.frame(stats_file = statsfiles,
                    FFT = qvalue1))

  message("Average FFT: ", mean(qvalue1))

}
