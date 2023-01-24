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
  id_request <- httr::GET("https://rest.uniprot.org/uniprotkb/stream",query=list(format="list",query=paste("model_organism:",taxonomy," AND reviewed:true AND keyword:KW-1185",sep="")))
  httr::stop_for_status(id_request,'download Uniprot ids from rest.uniprot.org')
  id_text <- httr::content(id_request,as='text')
  idlist <- unlist(strsplit(id_text,"\n"))
  data.env = getDataEnvironment()
  data.env[[ paste("gator.UniProtData.accs.",taxonomy,sep="") ]] <- idlist
  return (idlist)
}

#' @export
batchUpdateUniprotSequences <- function(taxonomy,reviewed=TRUE) {
  if (reviewed) {
    data_request <- httr::GET("https://rest.uniprot.org/uniprotkb/stream",query=list(format="tsv",fields="accession,sequence",query=paste("model_organism:",taxonomy," AND reviewed:true",sep="")),httr::progress())
  } else {
    data_request <- httr::GET("https://rest.uniprot.org/uniprotkb/stream",query=list(format="tsv",fields="accession,sequence",query=paste("model_organism:",taxonomy,sep="")),httr::progress())
  }
  httr::stop_for_status(data_request,'download Uniprot data from rest.uniprot.org')
  seq_data <- httr::content(data_request,as='text')
  target_file=tempfile()
  fileConn<-file(target_file)
  writeLines(seq_data, fileConn)
  close(fileConn)
  loadParsedJson('gator.UniProtData')
  target_env = getDataEnvironment()
  if (! exists('gator.UniProtData',target_env)) {
    target_env[['gator.UniProtData']] = data.frame( uniprot = character(0), sequence = character(0), stringsAsFactors=FALSE)
  }
  new_data = setNames(read.delim(target_file),c('uniprot','sequence'))
  new_data$uniprot = tolower(new_data$uniprot)
  assign('gator.UniProtData',rbind(get('gator.UniProtData',target_env),new_data),envir=target_env)
  writeParsedJson('gator.UniProtData')
  NULL
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
    return (unique(subset(gator.UniProtData, uniprot %in% tolower(accessions) )))
  }
  if (length(wanted_accs) > 100) {
    accgroups <- split(wanted_accs, ceiling(seq_along(wanted_accs)/100))
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

  accession <- toupper(paste(unlist(wanted_accs),collapse=","))

  fastas <- httr::GET("https://www.ebi.ac.uk/proteins/api/proteins",query=list(accession=accession), httr::add_headers(Accept = "text/x-fasta"))


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