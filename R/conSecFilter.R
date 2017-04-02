#' conSecFilter
#' @param peptideprofiles data.frame with protein_id, peptide_id + fraction number columns, filled with intensities
#' @param min_consecutive_ids minimum required consecutive ids to keep quantifications
#' @param remove_empty whether rows that are empty after filtering shall be removed (i.e. peptide entries with less than min_consecutive_ids)
#' @return filtered data.frame of same structure as input peptideprofiles table
conSecFilter<-function(peptideprofiles, min_consecutive_ids=3, remove_empty=TRUE)
  {
  # prepare input table
  peptide.traces <- as.matrix(peptideprofiles[,3:ncol(secdata)])
  rownames(peptide.traces) <- peptideprofiles$peptide_id
  # Define n as the number of columns
  labels <- peptideprofiles[, c("protein_id", "peptide_id")]
  data.filtered <- peptide.traces
  # Add 0-column
  data.filtered <- cbind(data.filtered, dummy=rep(0, nrow(data.filtered)))
  # Filter
  n <- ncol(data.filtered)
  # Set count variable to one
  tmp <- 1
  # Going through all rows, for all do:
  npeps <- nrow(data.filtered)
  for (x in 1:npeps) {
    message(paste('processed', x, 'of', npeps, 'rows'))
    tmp <- 1
    # Go through values in all cols, from left to right and ask the following:
    for (i in 1:n) {
      # Is value in row x, col i = 0?
      if (data[x,i] == 0) {
        # If the value is = o, tmp stays = 1
        tmp<-1
      }
      # If value is not = 0, then:
      else {
        # Look at next value at row x, col i+1, is it 0? If yes, then:
        if (data[x, i+1] == 0) {
          # Do nothing if the count variable is >3, i.e. there's 4 or more values in a stretch
          if (tmp > (min_consecutive_ids-1)) {}
          else {
            # Replace values if count variable is 3 or smaller at x, i, where the next value is 0. Replace with 0
            for (j in 0:tmp-1) {
              data[x, i-j] <- 0
            }
            tmp <- 1
          }
        }
        # If next value is not = 0, then count tmp+1
        else {
          tmp <- tmp+1
        }
      }
    }
  }

  # Remove dummy column
  data <- subset(data, select =-dummy)
  # Write data to traces object
  peptide.traces.filtered <- cbind(labels, data)
    if (remove_empty) {
    peptide.traces.filtered <- peptide.traces.filtered[rowSums(data) != 0, ]
  }
  return(peptide.traces.filtered)
}

