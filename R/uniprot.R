#' @importFrom httr GET
#' @importFrom httr content
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
  assign( paste("gator.UniProtData.accs.",taxonomy,sep=""), idlist, envir = .GlobalEnv )
  return (idlist)
}

#' @importFrom data.table rbindlist
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom plyr ldply
#' @importFrom data.table rbindlist
#' @export
getUniprotSequences <- function(accs,wait=0) {
  if (length(accs) > 200) {
    accgroups <- split(accs, ceiling(seq_along(accs)/200))
    accumulated_frame <- NULL
    pb <- txtProgressBar(min=0, max=length(accgroups),initial=0)
    for (i in seq_along(accgroups)) {
      frame <- suppressMessages(getUniprotSequences(as.vector(unlist(accgroups[i])),wait=5))
      if (is.null(accumulated_frame)) {
        accumulated_frame <- frame
      } else {
        accumulated_frame <- data.table::rbindlist(list(accumulated_frame,frame))
      }
      setTxtProgressBar(pb,i)
    }
    close(pb)
    return (unique(accumulated_frame))
  }
  wanted_accs <- accs
  cached <- loadParsedJson('UniProtData')
  if (dim(cached)[1] > 0) {
    assign("gator.UniProtData",cached, envir = .GlobalEnv)
  }
  cached <- NULL
  if (exists("gator.UniProtData")) {
    wanted_accs <- unique(wanted_accs[! wanted_accs %in% gator.UniProtData$uniprot ])
  } else {
    assign("gator.UniProtData",data.frame( uniprot = character(0), sequence = character(0), stringsAsFactors=FALSE), envir = .GlobalEnv)
  }
  if (length(wanted_accs) < 1) {
    return (subset(gator.UniProtData, uniprot %in% accs ))
  }
  message("Retrieving ",length(wanted_accs)," UniProt IDs")
  fastas <- httr::POST("http://www.uniprot.org/batch/",body=list(format='fasta',file=RCurl::fileUpload('upload',toupper(paste(unlist(wanted_accs),collapse="\n")))),multipart=TRUE)
  if (fastas$status_code != 200) {
    message("Could not retrieve ids")
    return(unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
  }
  contents <- httr::content(fastas,"text")
  seqs <- strsplit(sub("\n","\t", unlist(strsplit(contents,"\n>"))),"\t")
  seqs <- plyr::ldply(seqs,function(row) { c(  row[1] , gsub("\n","",row[2]) )  });
  seqs$V1 <- sub(">?[st][pr]\\|","",seqs$V1)
  seqs$V1 <- sub("\\|.*","",seqs$V1)
  if (length(names(seqs)) < 2) {
    message("Could not retrieve sequences")
    message(accs)
    return (unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
  }
  names(seqs) <- c('uniprot','sequence')
  seqs$uniprot <- tolower(seqs$uniprot)
  assign('gator.UniProtData', data.table::rbindlist( list(get('gator.UniProtData'), seqs) ), envir = .GlobalEnv)
  writeParsedJson(gator.UniProtData,'UniProtData')
  Sys.sleep(wait)
  return (unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
}