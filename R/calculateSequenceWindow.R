#'  Append a window column to a given data frame
#'
#'  This function silently drops any sites that do not have enough neighbouring
#'  amino acids to make a complete window.
#'
#'  @param  dataframe   data.frame with a uniprot column, and another numeric column with site
#'  @param  sitecol     Name of the column containing the site information
#'  @param  window      Half-width of window (\code{window} of \code{3} results in 7-aa windows)
#'  @examples
#'    df <- data.frame(uniprot='q14118',site=c(30,40,50))
#'    # 11-aa windows
#'    calculateSequenceWindow(df,'site',5)
#'    # 5-aa windows
#'    calculateSequenceWindow(df,'site',2)
#'  @export
calculateSequenceWindow <- function(dataframe,sitecol,window) {
  with_seqs <- unique(merge(dataframe,getUniprotSequences(unique(dataframe$uniprot)),by='uniprot'))
  with_seqs[[sitecol]] <- as.numeric(with_seqs[[sitecol]])
  with_seqs$endpos <- nchar(with_seqs$sequence) - with_seqs[[sitecol]] - window
  with_seqs$startpos <- with_seqs[[sitecol]] - window
  with_seqs <- subset(with_seqs, (endpos > 0 ))
  with_seqs <- subset(with_seqs, (startpos > 0 ))
  with_seqs$window <- substr(with_seqs$sequence,with_seqs[[sitecol]]-window,with_seqs[[sitecol]]+window)
  with_seqs$sequence <- NULL
  with_seqs$endpos <- NULL
  with_seqs$startpos <- NULL

  return(with_seqs)
}
