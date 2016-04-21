# @importFrom httr GET
# @importFrom httr content
#' Retrieve the set of UniProt accessions for a given taxonomy ID
#'
#' Only get the identifiers for entries where they are part of a complete proteome set (keyword 1185)
#' and are reviewed
#'
#' @param   taxonomy  NCBI taxonomy identifier to retrieve data for
#' @return  Vector of string accessions
#' @export
getUniprotIds <- function(taxonomy) {
  # We should add an option to get unreviewed ids here too
  if (exists(paste("gator.UniProtData.accs.",taxonomy,sep=""))) {
    return ( get(paste("gator.UniProtData.accs.",taxonomy,sep="")) )
  }
  message("Getting ID list")
  id_request <- httr::GET("http://www.uniprot.org/uniprot/",query=paste("query=taxonomy:",taxonomy,"+AND+reviewed:yes+AND+keyword:1185&force=yes&format=list",sep=""))
  id_text <- httr::content(id_request,as='text')
  idlist <- unlist(strsplit(id_text,"\n"))
  data.env = getDataEnvironment()
  data.env[[ paste("gator.UniProtData.accs.",taxonomy,sep="") ]] <- idlist
  return (idlist)
}

# @importFrom data.table rbindlist
# @importFrom httr POST
# @importFrom httr content
# @importFrom plyr ldply
# @importFrom data.table rbindlist
#' Retrieve UniProt sequences for a set of UniProt accessions
#' @param   accessions   Uniprot accessions to retrieve sequence data for
#' @return  Data frame with UniProt identifier and sequence
#' @export
getUniprotSequences <- function(accessions,wait=0) {
  wanted_accs <- gsub("-.*","",tolower(accessions))
  if (any(grepl("-",accessions))) {
    message("Isoform retrieval currently buggy in UniProt, removing isoforms")
  }
  loadParsedJson('gator.UniProtData')
  if (exists("gator.UniProtData")) {
    wanted_accs <- unique(wanted_accs[! wanted_accs %in% gator.UniProtData$uniprot ])
  } else {
    data.env = getDataEnvironment()
    data.env[[ 'gator.UniProtData']] <- data.frame( uniprot = character(0), sequence = character(0), stringsAsFactors=FALSE)
  }
  if (length(wanted_accs) < 1) {
    return (unique(subset(gator.UniProtData, uniprot %in% accessions )))
  }
  if (length(wanted_accs) > 200) {
    accgroups <- split(wanted_accs, ceiling(seq_along(wanted_accs)/200))
    accumulated_frame <- NULL
    pb <- txtProgressBar(min=0, max=length(accgroups),initial=0)
    for (i in seq_along(accgroups)) {
      frame <- suppressMessages(getUniprotSequences(as.vector(unlist(accgroups[i])),wait=5))
      if (is.null(accumulated_frame)) {
        accumulated_frame <- frame
      } else {
        accumulated_frame <- do.call(rbind,list(accumulated_frame,frame))
      }
      setTxtProgressBar(pb,i)
    }
    close(pb)
    return (unique(accumulated_frame))
  }

  message("Retrieving ",length(wanted_accs)," UniProt IDs")
  toupload <- tempfile()
  writeLines(toupper(paste(unlist(wanted_accs),collapse="\n")), toupload)
  acc_file <- httr::upload_file(toupload)

  fastas <- httr::POST("http://www.uniprot.org/batch/",body=list(format='fasta',file=acc_file),encode="multipart")

  unlink(toupload)

  if (fastas$status_code != 200) {
    message("Could not retrieve ids ",fastas$status_code)
    return(unique(subset(gator.UniProtData, uniprot %in% tolower(accessions) )))
  }
  contents <- httr::content(fastas,"text")
  seqs <- strsplit(sub("\n","\t", unlist(strsplit(contents,"\n>"))),"\t")
  seqs <- plyr::ldply(seqs,function(row) { c(  row[1] , gsub("\n","",row[2]) )  });
  seqs$V1 <- sub(">?[st][pr]\\|","",seqs$V1)
  seqs$V1 <- sub("\\|.*","",seqs$V1)
  if (length(names(seqs)) < 2) {
    message("Could not retrieve sequences")
    message(wanted_accs)
    return (unique(subset(gator.UniProtData, uniprot %in% tolower(accessions) )))
  }
  names(seqs) <- c('uniprot','sequence')
  seqs$uniprot <- tolower(seqs$uniprot)
  data.env = getDataEnvironment()
  data.env[[ 'gator.UniProtData']] <- as.data.frame(data.table::rbindlist( list(get('gator.UniProtData'), seqs)))
  gator.UniProtData$uniprot = tolower(gator.UniProtData$uniprot)
  writeParsedJson('gator.UniProtData')
  Sys.sleep(wait)
  return (unique(subset(gator.UniProtData, uniprot %in% tolower(accessions) )))
}