#' @export
calculateSequenceWindow <- function(dataframe,sitecol,window) {
  with_seqs <- merge(dataframe,gator.UniProtData,by='uniprot')
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
