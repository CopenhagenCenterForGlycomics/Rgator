# We need to do this to load up the libraries to start with
# library(devtools)
# load_all('~/dev/Rgator/')

#install_github("httr",username="hirenj")
#install_github("zeenogee/R-Websockets")


options(stringsAsFactors = FALSE)

gator.cache <- getOption("gator.alternativeCache")
gator.biocLite.stdlib <- FALSE

if (! is.null(gator.cache)) {
  dir.create(gator.cache,showWarnings=FALSE)
  if ( ! gator.cache %in% .libPaths() ) {
    .libPaths(c(.libPaths(), gator.cache ))
  }
} else {
  gator.cache <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(gator.cache,showWarnings=FALSE)
  gator.biocLite.stdlib <- TRUE
}

current_token <- NULL
connection_key <- NULL

#' @importFrom httr oauth_endpoint
#' @importFrom httr oauth_app
#' @importFrom httr oauth2.0_token
doSignin <- function() {
  if (!is.null(current_token)) {
    return((list(access_token=current_token)))
  }
  if (! is.null(getOption("current_token"))) {
    return(list(access_token=getOption("current_token")))
  }
	google <- httr::oauth_endpoint(NULL,"auth","token", base_url= "https://accounts.google.com/o/oauth2")
	gdrive <- httr::oauth_app("google",getOption("GoogleClientId"),secret=getOption("GoogleClientSecret"))

	token <- httr::oauth2.0_token(google, gdrive,scope="https://www.googleapis.com/auth/drive.readonly",use_oob = TRUE )

  current_token <- token$access_token
	return (token)
}

global_preferences <- NULL

getPreferences <- function(prefsId='0B5L9OYFFMK3dcmlBbUY3SXdHMk0') {
	if (is.null(global_preferences)) {
		return (getGoogleFile(prefsId))
	} else {
		return (global_preferences);
	}
}

#' @importFrom rjson fromJSON
#' @importFrom websockets websocket_write
askForSignin <- function(WS) {
  websockets::websocket_write(rjson::toJSON(list(message="upgradeConnection",data=getOption("connection_key"))), WS)
}

acceptToken <- function(json) {
  message("Successfully connected to GlycoDomain viewer")
  options(current_token = json$data$authtoken)
  options(connection_key = json$data$connectionkey)
}

#' @importFrom websockets daemonize
#' @importFrom websockets create_server
#' @importFrom websockets websocket_write
#' @importFrom websockets websocket_close
#' @importFrom rjson fromJSON
#' @importFrom websockets setCallback
#' @export
gatorConnector <- function() {
  server = websockets::create_server(port=8880)
  websockets::daemonize(server)

  view <- function(prot) {
    if (length(server$client_sockets) > 0) {
      websockets::websocket_write(rjson::toJSON(list(message="showProtein", data=unique(prot))), server$client_sockets[[1]])
    }
  }
  align <- function(prot) {
    if (length(server$client_sockets) > 0) {
      websockets::websocket_write(rjson::toJSON(list(message="alignProtein", data=unique(prot))), server$client_sockets[[1]])
    }
  }
  compactRenderer <- function() {
    if (length(server$client_sockets) > 0) {
      websockets::websocket_write(rjson::toJSON(list(message="compactRenderer", data="")), server$client_sockets[[1]])
    }
  }

  receiver <- function(DATA,WS,...) {
    json <- rjson::fromJSON(rawToChar(DATA))
    if (json$message == "token") {
      acceptToken(json)
    }
  }

  websockets::setCallback("receive",receiver,server)
  websockets::setCallback("established", askForSignin, server)

  stopConnector <- function() {
    websockets::websocket_close(server)
    options(current_token = NULL)
    options(connection_key = NULL)
    rm("stopConnector",envir=.GlobalEnv)
    rm("Viewp",envir=.GlobalEnv)
    rm("alignProtein",envir=.GlobalEnv)
    rm("compactRenderer",envir=.GlobalEnv)
  }
  assign("Viewp",view,envir = .GlobalEnv)
  assign("alignProtein",align,envir = .GlobalEnv)
  assign("compactRenderer",compactRenderer,envir=.GlobalEnv)
  assign("stopConnector",stopConnector,envir=.GlobalEnv)
}

#' @importFrom rjson fromJSON
#' @export
syncDatasets <- function() {
  files <- list.files(pattern = "\\.domaintoolsession$")
  all_prefs <- list()
  if (length(files) > 0) {
    all_prefs <- lapply(files,function(filename) {
      fileConn <- file(filename,"r")
      data <- rjson::fromJSON(,fileConn)
      close(fileConn)
      return (data)
    })
  } else {
     all_prefs = c(getPreferences('0By48KKDu9leCeXN3cEhYZGlwVjQ'))
  }

  lapply(all_prefs,function(prefs) {
    lapply(seq_along(prefs$user_datasets), function(i) { key <- names(prefs$user_datasets)[i]; config <- prefs$user_datasets[[key]]; downloadDataset(key,config) } )
  });

}

#' @importFrom plyr ldply
jsonParser <- function(data,keys) {
  currkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key))>0) { strsplit(key,".",fixed=TRUE)[[1]][1] } }))
  currkeys <- unlist(currkeys[!sapply(currkeys,is.null)])

  localkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key)) < 1) { return(key) } }))
  localkeys <- unlist(localkeys[!sapply(localkeys,is.null)])

  keys <- unlist(lapply(keys,FUN=function(key) { sub("[^\\.]+\\.","",key) }))
  if (length(currkeys) == 1 && currkeys[1] != "*" && (length(keys) > length(currkeys) || keys[1] != currkeys[1]) ) {
    kidframe <- (plyr::ldply(data[[ currkeys[1] ]],.fun=function(dat) { parsed <- jsonParser(dat,keys[!keys %in% localkeys]); return(parsed); }))
    ## Sometimes if the sites array was empty, we were not repeating the other fields
    ## { "alwayspresent" : "ABCD", "sites" : [] } , { "alwayspresent": "ABCDE" , "sites" : [ { "somekey" : "ABC" , "someotherkey" : "DEF"}]}
    if (dim(kidframe)[2] < length(keys[!keys %in% localkeys])) {
      kidframe <- data.frame(lapply( keys[!keys %in% localkeys],function(el) { return (c(NA)) }))
      names(kidframe) <- c(1: length( keys[!keys %in% localkeys] ) )
    }
    for (local in localkeys) {
      if (length(data[[local]]) < 2) {
        kidframe[[local]] <- rep(data[[local]],dim(kidframe)[1])
      } else {
        subframe <- data.frame(lapply( data[[local]] , function(val) {
          rep(val,dim(kidframe)[1])
        }))
        names(subframe) <- c(1:length(data[[local]]))
        kidframe <- cbind(kidframe,subframe)
      }
    }
    return (kidframe)
  } else {
    if (length(keys) == 1) {
      if (length(currkeys) > 0 && currkeys[1] == "*") {
        retval <-lapply( data, FUN=function(dat) { return ((jsonParser(dat,keys))) }  )
        return (plyr::ldply(retval))
      }
      ##
      ## { "orthologs" : {"TAX1" : ["a","b","c"], "TAX2" : ["a","b"] }}
      ##
      if (length(unique(lapply(data[[keys[1]]],length))) > 1) {
        return (plyr::ldply(data[[keys[1]]],.fun=function(list) { return(data.frame(V1=list)) }))
      }
      return (plyr::ldply(data[[keys[1]]]))
    }
    datlist <- sapply(c(1:(length(keys)-1)),FUN=function(keyidx) {
      key <- keys[keyidx]
      if (is.null(data[[key]]) || length(data[[key]]) == 0) {
        data[[key]] = NA
      }
      return (data[[key]])
    })
    key <- keys[length(keys)]
    if (is.null(data[[key]]) || length(data[[key]]) == 0) {
      data[[key]] = c(NA)
    }
    datlist <- c( datlist, list(  data[[ keys[length(keys)]  ]] ))

    retval <- do.call("cbind",datlist)
    return (retval)
  }
}

#' @export
downloadOrthologies <- function() {
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/homologene/',list(type='gatorURL',title='orthology.homologene'))
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/orthology.treefam.9/',list(type='gatorURL',title='orthology.treefam'))
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/orthology.inparanoid.8_1/',list(type='gatorURL',title='orthology.inparanoid'))
  assign('gator.orthology',rbind(gator.homologene,`gator.orthology.treefam`,`gator.orthology.inparanoid`),envir=.GlobalEnv)
}

#' @export
downloadDomains <- function(organism=c('9060')) {
  if ("9060" %in% organism) {
    downloadDataset('http://glycodomain-data.glycocode.com/data/latest/fulldomains/',list(type='gatorURL',title='fulldomains'))
  }
  downloadDataset(paste('http://glycodomain-data.glycocode.com/data/latest/domains.',organism,'/',sep=''),list(type='gatorURL',title=paste('domains.',organism,sep='')))
}

#' @importFrom plyr ldply
#' @importFrom data.table rbindlist
downloadDataset <- function(set,config,accs=c(),etagcheck=TRUE) {
  message("Downloading ",set)
  if (length(grep("gatorURL",config[['type']])) > 0) {
    # This is only the key for the hash in the prefs file
    # We need to be storing the type with the key
    # set.type == gatorURL
    # If there's no trailing slash - we want to do a query
    if (substr(set, nchar(set), nchar(set)) != '/' || length(accs) > 0) {
      if (length(accs) > 50) {
        accgroups <- split(accs, ceiling(seq_along(accs)/20))
        accumulated_frame <- NULL
        for (i in seq_along(accgroups)) {
          message("Retrieving for ",set," step ",i," out of ",length(accgroups))
          frame <- downloadDataset(set,config,accs=accgroups[i],etagcheck=FALSE)
          if (is.null(accumulated_frame)) {
            accumulated_frame <- frame
          } else {
            accumulated_frame <- data.table::rbindlist(list(accumulated_frame,frame))
          }
          if ( i == length(accgroups)) {
            assign(paste("gator.",gsub("[[:space:]]|-","_",config[['title']]),sep=""),accumulated_frame, envir = .GlobalEnv)
          }
        }
        return (frame)
      } else {
        data <- getGatorSnapshotSubset(set,accs)
      }
    } else {

      # Otherwise we should switch to grabbing a subset
      data <- getGatorSnapshot(substr(set, 1, nchar(set)-1),gsub("[[:space:]]|-","_",config[['title']]))
    }
  } else {
    if (config[['type']] == 'dataset') {
      # set.type == dataset
      data <- getGoogleFile(set)
    } else {
      message("Don't know how to retrieve ",set)
      return()
    }
  }

  # We should check to make sure that the etag for the current data frame in the global namespace
  # matches with the etag for the current data, so we can find out if we have to redo the parsing of
  # the data
  message("Checking etags of pre-parsed data")
  if(etagcheck) {
    if (exists( paste("gator.",gsub("[[:space:]]|-","_",data$title),sep="") )) {
      message("We have data loaded up in the environment")
      frame <- get( paste("gator.",gsub("[[:space:]]|-","_",data$title),sep="") )
      if (!is.null(attributes(frame)$etag) && attributes(frame)$etag == format(data$etag,scientific=FALSE)) {
        return ()
      }
      message("No etag match, continuing")
    } else {
      frame <- loadParsedJson(data$title)
      if (!is.null(attributes(frame)$etag) && attributes(frame)$etag == format(data$etag,scientific=FALSE)) {
        message("We have data that has already been parsed")
        assign(paste("gator.",gsub("[[:space:]]|-","_",data$title),sep=""),frame, envir = .GlobalEnv)
        return ()
      }

    }
  }

  message("Preparing to parse data")
  if (!is.null(data$defaults$rKeys) && data$defaults$rKeys == c('base64')) {
    all_prots <- names(data$data)
    names(all_prots) <- all_prots
    frame <- llply(all_prots,.fun=function(uprot) {
      if (nchar(data$data[[uprot]]) > 0) {
        decodeBase64(data$data[[uprot]])
      }
    },.progress="text")
    assign(paste("gator.",gsub("[[:space:]]|-","_",data$title),sep=""),frame, envir = .GlobalEnv)
    attributes(frame)$etag <- data$etag
    writeParsedJson(frame,data$title)
    return (data.frame())
  } else if (!is.null(data$defaults$rKeys) && length(data$defaults$rKeys) > 0) {
    all_prots <- names(data$data)
    frame <- data.table::rbindlist(llply(all_prots,.fun=function(uprot) {
      frm <- jsonParser(data$data[[uprot]],data$defaults$rKeys )

      # We should get a data frame out from the jsonParser - attach the uniprot id as
      # another column into the data frame

      frm$uniprot <- rep(uprot,dim(frm)[1])
      return(frm)
    },.progress="text"))

    # We need to re-arrange the columns here so that the uniprot column
    # ends up as the first column for consistency
    wanted_cols <- names(frame)
    frame <- subset(frame,select=c('uniprot',wanted_cols[!wanted_cols == 'uniprot']))
    setnames(frame, c('uniprot', data$defaults$rKeys, rep('NA',dim(frame)[2] - (length(data$defaults$rKeys)+1))))
  } else {
    # Assume that we're just pulling out sites from the data sets if we're not given a particular
    # key to iterate over
    frame <- plyr::ldply(data$data,.fun=function(dat) { return (plyr::ldply(dat$sites)) },.progress="text")
    currnames <- names(frame)
    currnames[1] <- 'uniprot'
    names(frame)<- currnames

  }
  if (!is.null(data$defaults$rNames)) {
    names(frame) <- c('uniprot',data$defaults$rNames)
  }
  attributes(frame)$etag <- data$etag

  assign(paste("gator.",gsub("[[:space:]]|-","_",data$title),sep=""),frame, envir = .GlobalEnv)

  writeParsedJson(frame,data$title)

  frame
}

loadParsedJson <- function(title) {
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  frame <- data.frame()

  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    origData <- load(fileConn)
    close(fileConn)
  }
  return(frame)
}

writeParsedJson <- function(frame,title) {
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  save(frame,file=filename)
}

#' @importFrom RCurl base64Decode
decodeBase64 <- function(base64) {
  decoded <- RCurl::base64Decode(base64,'raw');
  readBin( decoded, 'numeric', length(decoded)/4 ,4);
}

#' @importFrom rjson fromJSON
testParseBJson <- function(filename) {
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    message("Loading cached file")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  all_prots <- names(origData$data)
  names(all_prots) <- all_prots
  frame <- llply(all_prots,.fun=function(uprot) {
    if (nchar(origData$data[[uprot]]) > 0) {
      decodeBase64(origData$data[[uprot]])
    }
  },.progress="text")

  return (frame)
}

#' @importFrom rjson fromJSON
#' @importFrom data.table rbindlist
testParseJson <- function(filename) {
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    message("Loading cached file")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  all_prots <- names(origData$data)
  frame <- data.table::rbindlist(llply(all_prots,.fun=function(uprot) {
    frm <- jsonParser(origData$data[[uprot]],origData$defaults$rKeys )

    # We should get a data frame out from the jsonParser - attach the uniprot id as
    # another column into the data frame

    frm$uniprot <- rep(uprot,dim(frm)[1])
    return(frm)
  },.progress="text"))

  # We need to re-arrange the columns here so that the uniprot column
  # ends up as the first column for consistency

  wanted_cols <- names(frame)
  frame <- subset(frame,select=c('uniprot',wanted_cols[!wanted_cols == 'uniprot']))
  setnames(frame, c('uniprot', origData$defaults$rKeys, rep('NA',dim(frame)[2] - (length(origData$defaults$rKeys)+1))))

  if (!is.null(origData$defaults$rNames)) {
    names(frame) <- c('uniprot',origData$defaults$rNames)
  }

  return (frame)
}

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

calculatePWM <- function(dataframe,windowcol,codes=c('A','C', 'D','E','F','G','H','I','J','K','L','M','N','P','Q','R','S','T','V','W','Y','Z')) {
  aas <- t(data.matrix(apply(as.array(dataframe[[windowcol]]),1,FUN=function(x) { unlist(strsplit(x,''))})))
  freqs <- apply(aas,2,function(x) { table(x) })
  sapply(1:dim(aas)[2],function(pos) {
    pos_freqs <- freqs[[pos]];
    total <- sum(pos_freqs)
    sapply(codes,function(aa) {
      if (! aa %in% names(pos_freqs) ) {
        return (0)
      }
      val <- pos_freqs[names(pos_freqs) == aa]
      if ( ! is.null(val) ) {
        return(as.integer(val) / total)
      }
    });
  })
}

#' @importFrom plyr .
#' @importFrom plyr ddply
#' @export
getDomainSets <- function( inputsites, sitecol, domaindata, max_dom_proportion=0.81, stem_distance=100  ) {
  message("Retrieving Uniprot sequences")
  seqdat <- getUniprotSequences(unique(inputsites$uniprot))
  message("Retrieved sequences")
  domdat <- merge(merge(domaindata,subset(inputsites,! is.na(inputsites[[sitecol]])),by='uniprot',allow.cartesian=TRUE),seqdat,by='uniprot',allow.cartesian=TRUE)
  domdat$aalength <- nchar(domdat$sequence)
  domdat$sitekey <- paste(domdat$uniprot,'-',domdat[[sitecol]],sep='')
  domdat$sequence <- NULL
  real <- subset( domdat, ((as.numeric(end) - as.numeric(start)) / aalength < max_dom_proportion ))
  inside <- subset( subset (  real ,  ( (as.numeric(real[[sitecol]]) >= as.numeric(start)) & (as.numeric(real[[sitecol]]) <= as.numeric(end))  )  ), ! grepl("tmhmm",dom))
  outside <- subset ( real , ! sitekey %in% inside$sitekey & dom != 'tmhmm-outside' & dom != 'tmhmm-inside' )
  sitekeys_nterm <- unique(subset( outside, as.numeric(outside[[sitecol]]) < as.numeric(start) )$sitekey)
  sitekeys_cterm <- unique(subset( outside, as.numeric(outside[[sitecol]]) > as.numeric(end) )$sitekey)
  between <- subset( outside, sitekey %in% intersect(sitekeys_nterm, sitekeys_cterm ) )
  message("Identifying type II proteins")
  typeii <- unique(plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- unique(subset(input,select=c('dom','start','end','sitekey')))
    signals <- subset(df, dom == 'SIGNALP')
    tms <- subset(df, dom == "tmhmm-TMhelix")
    if (dim(signals)[1] < 1) {
      if (dim(tms)[1] == 1) {
        return (input)
      }
      return ()
    }
    signal_start <- signals$start[1]
    signal_end <- signals$end[1]
    tms <- subset(df, dom == "tmhmm-TMhelix" & as.numeric(start) < as.numeric(signal_end) )
    other_tms <- subset(df, dom == "tmhmm-TMhelix" & as.numeric(start) > as.numeric(signal_end) )

    if (dim(tms)[1] > 0 & dim(other_tms)[1] == 0) {
      return (input)
    }
    return ()
  },.progress="text")$uniprot)
  message("Identifying type II stems")
  stem_typeii <- plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the N-terminal domains, our site is C-terminal of the domain
    filtered <- subset(df,siteend>0)
    filtered <- filtered[order(filtered$siteend),]
    # We have a SIGNALP or TMHELIX -- something...
    # This is a type II transmembrane or it is secreted
    if (( (filtered$dom[1] == "SIGNALP") | grepl("tmhmm-TMhelix",filtered$dom[1]) ) & (filtered$uniprot[1] %in% typeii )) {
      wanted <- subset(df, ((startsite > 0 & startsite <= stem_distance) | (siteend > 0 & siteend <= stem_distance & siteend <= filtered$siteend[1] )))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
  },.progress="text")
  message("Identifying Signalp stems")
  signalp_stem <- plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the N-terminal domains, our site is C-terminal of the domain
    filtered <- subset(subset(df,siteend>0),dom != 'tmhmm-TMhelix')
    filtered <- filtered[order(filtered$siteend),]
    # We have a SIGNALP -- something...
    # If we've got a signalp then it is secreted


    if ( dim(filtered)[1] > 0 & ( (filtered$dom[1] == "SIGNALP") ) & (! filtered$uniprot[1] %in% typeii )) {
      # Anything with a transmembrane helix is not secreted
      wanted <- subset(df, ((siteend > 0 & siteend <= stem_distance & siteend <= filtered$siteend[1] )))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
  },.progress="text")
  message("Identifying type I stems")
  stem_typei <- plyr::ddply(subset(between, ! sitekey %in% signalp_stem$sitekey),plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the C-terminal domains, our site is N-terminal of the domain
    # This is a type I transmembrane
    filtered <- subset(df,startsite>0)
    filtered <- filtered[order(filtered$startsite),]
    # We have a something -- TMHELIX

    # If we have a domain spanning the transmembrane, we want to call this a type I transmembrane too

    if ( ((dim(filtered)[1] >= 1) & ( filtered$dom[1] == "tmhmm-TMhelix" )) | ((dim(filtered)[1] >= 2) & (filtered$dom[2] == "tmhmm-TMhelix") & (as.numeric(filtered$end[1]) > as.numeric(filtered$end[2]) ) )) {
      wanted <- subset(df, ((startsite > 0 & startsite <= stem_distance & startsite <= filtered$startsite[1] ) | (siteend > 0 & siteend <= stem_distance)))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
    return ()
  },.progress="text")

  interdomain <- subset(between, ! sitekey %in% stem_typei$sitekey & ! sitekey %in% stem_typeii$sitekey & ! sitekey %in% signalp_stem$sitekey )
  norc <- subset(outside, ! sitekey %in% between$sitekey )
  #Stem = Betweeen where closest N-terminal = SIGNALP/TMHMM
  #ddply between by sitekey if (site - end), sort asc [1] $dom == tmhmmm/signalp return df
  #                         if (start - site), sort asc [1] $dom == tmhmm/signalp return df
  #                         else return empty
  return ( list( all=domdat, real=real, inside=inside, outside=outside, between=between, stem=rbind(stem_typei,stem_typeii,signalp_stem), stem.typeii=stem_typeii, stem.typei=stem_typei, stem.signalp=signalp_stem, interdomain=interdomain, norc=norc  )  )
}

cacheFile <- function(url,fileId,gzip=F) {
  filename <- file.path(gator.cache,paste("gator-",fileId,sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    if (gzip) {
      filename <- gzfile(filename)
    }
    return (read.delim(filename,header=F,sep='\t'))
  }
  download.file(url,filename)

  if (gzip) {
    filename <- gzfile(filename)
  }
  return (read.delim(filename,header=F,sep='\t'))
}

#' @export
cddidToSuperfamily <- function(cddids) {
  cdd_superfamily_links <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links","cdd-family-superfamily")[,c(1,3)]
  names(cdd_superfamily_links) <- c('cddid','clusterid')
  for(i in 1:nrow(cdd_superfamily_links)) {
    cddids[ cddids==cdd_superfamily_links[[i,'cddid' ]] ] <- cdd_superfamily_links[[ i , 'clusterid'  ]]
  }
  cddids
}

#' @export
getCddNames <- function(cddids) {
  cddid_all <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz","cddid-all",gzip=T)
  names(cddid_all) <- c('id','dom','short','long')
  merge(data.frame(dom=cddids),cddid_all,by='dom',all.x=T)
}

uniqueframe <- function(set){
  return(as.data.frame(unique(as.matrix(set))))
}

#' @export
generateLogoPlot <- function(dataframe,windowcol,frequencies=c()) {
  uniprot_2013_12_freq <- list(A=0.0825,R=0.0553,N=0.0406,D=0.0545,C=0.0137,Q=0.0393,E=0.0675,G=0.0707,H=0.0227,I=0.0595,L=0.0966,K=0.0584,M=0.0242,F=0.0386,P=0.0470,S=0.0657,T=0.0534,W=0.0108,Y=0.0292,V=0.0686)
  if(length(frequencies) < 1) {
    frequencies <- uniprot_2013_12_freq
  }
  pwm <- calculatePWM(dataframe,windowcol,names(frequencies))
  return(berrylogo(pwm,frequencies))
}

#' @importFrom grid gTree
#' @importFrom grid gList
#' @importFrom grid grid.grabExpr
#' @importFrom grid grid.draw
#' @importFrom VennDiagram venn.diagram
#' @export
generateVennDiagram <- function(data=list(),title="Venn Diagram") {
  require(VennDiagram)
  return (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(VennDiagram::venn.diagram(data,filename=NULL,main=title)))), cl=c("arrange", "ggplot")))
}

#' @importFrom plyr laply
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 coord_trans
#' @importFrom ggplot2 theme_bw
berrylogo<-function(pwm,backFreq,zero=.0001){
  pwm[pwm==0]<-zero
  bval <- plyr::laply(names(backFreq),function(x) {  pwm[x,] / backFreq[[x]] })
  row.names(bval)<-names(backFreq)
  hydrophobicity <- c(rep('#0000ff',6),rep('#00ff66',6),rep('#000000',8))
  names(hydrophobicity) <-  unlist(strsplit('RKDENQSGHTAPYVMCLFIW',''))
  chemistry <- c(rep('#00ff66',7),rep('#0000ff',3),rep('#ff0000',2),rep('#000000',8))
  names(chemistry) <- unlist(strsplit('GSTYCQNKRHDEAVLIPWFM',''))
  bval <- bval[names(hydrophobicity),]
  window_size = floor( 0.5*length(dimnames(bval)[[2]]) )
  dimnames(bval)[[2]]<- c((-1*window_size):window_size)
  p<-ggplot2::ggplot(reshape2::melt(bval,varnames=c("aa","pos")),ggplot2::aes(x=pos,y=value,label=aa))+
    ggplot2::geom_line(ggplot2::aes(y=1), colour = "grey",size=2)+
    ggplot2::geom_text(ggplot2::aes(colour=factor(aa)),face='bold',size=8)+scale_colour_manual(values=chemistry)+
    ggplot2::theme(legend.position="none")+
    ggplot2::scale_x_continuous(name="Position",breaks=(-1*window_size):window_size)+
    ggplot2::scale_y_continuous(name="Relative frequency",breaks=c(0.0625,0.125,0.25,0.5,1,seq(2,20,by=2)),limits=c(2^-11,20),label=function(x) format(x,nsmall = 2,drop0trailing=T,scientific = FALSE))+
    ggplot2::coord_trans(y='log2')+
    ggplot2::theme_bw()
  return(p)
}

#' Get entrez gene identifiers for a set of UniProt ids
#'
#' @importFrom AnnotationDbi toTable
#' @export
getEntrezIds <- function(organism,ids) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  uprotmap <- sub("\\.db","UNIPROT",dbname)
  uprots <- unique(intersect(  toupper(ids),mappedRkeys(get(uprotmap,asNamespace(dbname)))))
  entrez_ids <- AnnotationDbi::toTable(revmap(get(uprotmap,asNamespace(dbname)))[  uprots ])
  if ("systematic_name" %in% names(entrez_ids)) {
    entrez_ids <- subset(merge(entrez_ids,AnnotationDbi::toTable(get(sub("\\.db","ENTREZID",dbname),asNamespace(dbname))[]),by='systematic_name'),select=c('gene_id','uniprot_id'))
  }
  names(entrez_ids) <- c('gene_id','uniprot')
  gene_ids <- entrez_ids$gene_id
  names(gene_ids) <- entrez_ids$uniprot
  return (gene_ids)
}

#' Get Gene names for a set of UniProt identifiers
#'
#' @export
getGeneNames <- function(organism,ids) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  library(dbname,character.only=TRUE)
  wanted_cols <- c('UNIPROT', 'SYMBOL' )
  if (as.character(organism) == '4932') {
    wanted_cols <- c('UNIPROT','GENENAME')
  }
  mapping <- select(get(dbname,asNamespace(dbname)), keys=c(toupper(ids)), columns=wanted_cols, keytype="UNIPROT")
  names(mapping) <- c('uniprot','symbol')
  return (mapping)
}


#' Convert entrez ids to a data frame of uniprot, genename and entrez id
#'
#' @importFrom AnnotationDbi select
#' @export
convertEntrezIds <- function(organism,ids=c()) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  library(dbname,character.only=TRUE)
  if (length(ids) < 1) {
    return( data.frame(geneid=NA,uniprot=NA,genename=NA)[numeric(0),] )
  }
  wanted_cols <- c('UNIPROT', 'SYMBOL', 'ENTREZID')
  if (as.character(organism) == '4932') {
    wanted_cols <- c('UNIPROT', 'GENENAME', 'ENTREZID')
  }
  retdata <- AnnotationDbi::select(get(dbname,asNamespace(dbname)), keys=as.character(ids), columns=wanted_cols, keytype="ENTREZID")
  names(retdata) <- c('geneid','uniprot','genename')
  retdata
}

#' @export
getGOenrichmentGenes <- function(enrichment,wanted_terms=c(),organism=9606) {
  unique(convertEntrezIds(organism,unique(unlist(sapply(wanted_terms, function(go) {   geneIdsByCategory(enrichment)[[go]]  })))))
}

#' Get the GO terms associated with the given UniProt ids
#'
#' @importFrom GO.db GO.db
#' @importFrom AnnotationDbi toTable
#' @export
getGOTerms <- function(organism,uniprots) {
  getBiocLiteLib('GO.db')
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  query_ids <- getEntrezIds(organism,uniprots)
  if (as.character(organism) == '4932') {
    entrez_ids <- query_ids
    entrez_mapping <- AnnotationDbi::toTable(revmap(org.Sc.sgd.db::org.Sc.sgdENTREZID)[query_ids])
    query_ids <- unlist(entrez_mapping[1])
  }
  godb <- get( sub("\\.db","GO",organisms[as.character(organism)] ), envir=get(unlist( organisms[as.character(organism)])))
  terms <- AnnotationDbi::toTable(godb[query_ids])
  names(terms) <- c('gene_id','go_id','Evidence','Ontology')
  if (as.character(organism) == '4932') {
    names(terms) <- c('systematic_name','go_id','Evidence','Ontology')
    terms <- merge(terms, data.frame(systematic_name=query_ids), by='systematic_name')
    terms <- merge(terms, entrez_mapping, by='systematic_name')
    terms <- uniqueframe( merge(terms, data.frame(uniprot=names(entrez_ids),gene_id=entrez_ids) , by='gene_id'))
    terms$uniprot <- tolower(terms$uniprot)
  } else {
    terms <- merge( terms, data.frame(gene_id=query_ids,uniprot=tolower(names(query_ids))), by='gene_id')
  }
  if (dim(terms)[1] > 0) {
    go_terms <- (select(GO.db::GO.db, terms$go_id,"TERM"))
    names(go_terms) <- c('go_id','term')
    terms <- merge( terms, go_terms, by='go_id')
    terms <- uniqueframe(terms)
    return (terms)
  }
  return (terms)
}

#' @importFrom AnnotationDbi toTable
#' @export
GO_parents <- function(node = "GO:0008150", ontology = "BP") {
    if (ontology == "BP") GOPARENTS <- GO.db::GOBPPARENTS
    if (ontology == "CC") GOPARENTS <- GO.db::GOCCPARENTS
    if (ontology == "MF") GOPARENTS <- GO.db::GOMFPARENTS
    kids <- node

    # initialize output
    out <- c(kids)
    # do the following until there are no more parents
    while (any(!is.na(kids))) {
        # Get the unique children of the parents (that aren't NA)
        #children <- unique(unlist(mget(parents[!is.na(parents)], envir=GOCHILDREN)))
        parents <- unique(unlist(AnnotationDbi::toTable(GOPARENTS[kids[!is.na(kids)]])[2]))
        # append chldren to beginning of `out`
        # unique will keep the first instance of a duplicate
        # (i.e. the most recent child is kept)
        out <- unique(append(parents[!is.na(parents)], out))

        # children become the parents of the next generation
        kids <- parents
        if ('all' %in% kids) {
          kids <- c(NA)
        }
    }
    return(out)
}

#' @importFrom AnnotationDbi toTable
#' @export
GO_children <- function(node = "GO:0008150", ontology = "BP") {
    if (ontology == "BP") GOCHILDREN <- GO.db::GOBPCHILDREN
    if (ontology == "CC") GOCHILDREN <- GO.db::GOCCCHILDREN
    if (ontology == "MF") GOCHILDREN <- GO.db::GOMFCHILDREN
    parents <- node

    # initialize output
    out <- c(parents)
    # do the following until there are no more parents
    while (any(!is.na(parents))) {
        # Get the unique children of the parents (that aren't NA)
        #children <- unique(unlist(mget(parents[!is.na(parents)], envir=GOCHILDREN)))
        children <- unique(unlist(AnnotationDbi::toTable(GOCHILDREN[parents[!is.na(parents)]])[1]))

        # append chldren to beginning of `out`
        # unique will keep the first instance of a duplicate
        # (i.e. the most recent child is kept)
        out <- unique(append(children[!is.na(children)], out))

        # children become the parents of the next generation
        parents <- children
    }
    return(out)
}

getBiocLiteLib <- function(dbname) {
  if ( ! library(dbname,character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    if ( gator.biocLite.stdlib) {
      biocLite(dbname)
      return
    }
  }

  if ( ! suppressMessages(library(dbname,lib.loc=c(gator.cache),character.only=TRUE,logical.return=TRUE,quietly=TRUE)))  {
    if (! gator.biocLite.stdlib) {
      biocLite(dbname,lib=gator.cache)
      return
    }
  }

}

#' @importFrom AnnotationDbi toTable
#' @export
getGOEnrichment <- function(organism,uniprots,query_ids=c(),universe=c(),ontology='BP',direction='over',supplemental.terms=NA,conditional=TRUE) {
  getBiocLiteLib("GO.db")
  getBiocLiteLib("GOstats")
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  library(dbname,character.only=TRUE)
  if (length(query_ids) < 1 & length(uniprots) > 0) {
    query_ids <- getEntrezIds(organism,uniprots)
  }
  if (as.character(organism) == '4932' & length(query_ids) > 0) {
    entrez_ids <- query_ids
    entrez_mapping <- AnnotationDbi::toTable(revmap(org.Sc.sgdENTREZID)[query_ids])
    query_ids <- unlist(entrez_mapping[1])
  }

  if (length(universe) < 1) {
    universe <- as.list( mappedkeys(get( sub("\\.db","GO",dbname ),asNamespace(dbname) )) )
  } else {
    universe <- getEntrezIds(organism,universe)
  }
  library('GOstats')
  if (is.na(supplemental.terms)) {
  params <- new('GOHyperGParams',
              geneIds=query_ids,
              universeGeneIds=unlist(universe),
              ontology=ontology,
              pvalueCutoff=0.05,
              conditional=conditional,
              testDirection=direction,
              annotation=as.character(organisms[as.character(organism)])
             )
  } else {
    library("GSEABase")
    super_terms = subset( unique(rbind ( as.data.frame(AnnotationDbi::toTable( get( sub("\\.db","GO", dbname),asNamespace(dbname) )  )), supplemental.terms)), Ontology==ontology)
    frame=GOFrame(data.frame(super_terms$go_id, super_terms$Evidence, super_terms$gene_id),organism=as.character(organism))
    allFrame=GOAllFrame(frame)
    gsc <- GeneSetCollection(allFrame, setType = GOCollection())
    params <- GSEAGOHyperGParams(name="Supplemental GO annotation",
                             geneSetCollection=gsc, geneIds = query_ids, universeGeneIds = unlist(universe),
                             ontology = ontology, pvalueCutoff = 0.05, conditional = conditional, testDirection =  direction)
  }
  hgOver <- hyperGTest(params)
  return (hgOver)
}

#' Get up to top level GO terms
#' @importFrom data.table rbindlist
#' @importFrom AnnotationDbi Term
#' @export
tableGOTerms <- function(organism,uniprot,wanted=c('GO:0030054','GO:0005886','GO:0012505','GO:0005783','GO:0005794','GO:0005764','GO:0016020','GO:0005739','GO:0005634','GO:0005576','GO:0005773','GO:0005829'),ontology='CC') {
  all_terms <- getGOTerms(organism,uniprot)
  go_ids <- sapply( wanted, function(x) { GO_children(x,ontology) }, USE.NAMES=T,simplify=F )
  data.table::rbindlist(sapply(names(go_ids),function(x) { expand.grid(go_id=x,term=AnnotationDbi::Term(x),uniprot=unique(subset(all_terms, go_id %in% go_ids[[x]] & uniprot %in% uniprot )$uniprot)) },simplify=F))
}

#' @importFrom plyr ddply
#' @importFrom plyr .
addSiteColumn <- function(dataframe) {
  dataset <- plyr::ddply(dataframe,plyr::.(uniprot),function(df) { vals <- c(1:(nrow(df))); df$site <- vals; return (df); })
  eval.parent(substitute(dataframe<-dataset))
}

#' @importFrom rjson fromJSON
#' @importFrom httr GET
#' @importFrom httr content
getGatorSnapshotSubset <- function(fileId,accs) {
  url <- paste('http://localhost:3001/data/history/',fileId,'?accs=',paste(tolower(unlist(accs)),collapse=','),sep='')
  config <- list()
  file_request <- httr::GET(url,config=config)
  if (file_request$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return ()
  }
  if (file_request$status_code > 400 && file_request$status_code < 500) {
    message("Could not retrieved data from: ",url," got status code ",file_request$status_code)
    return ()
  }

  message(file_request$status_code)
  message("Retrieving data from Gator for ",httr::content(file_request)$title)
  retval <- rjson::fromJSON(httr::content(file_request,"text"))
  retval$etag <- format(retval$etag,scientific=FALSE)
  return (retval)
}

#' @importFrom rjson fromJSON
#' @importFrom httr GET
#' @importFrom httr content
getGatorSnapshot <- function(gatorURL,fileId) {
  filename <- file.path(gator.cache,paste("gator-",fileId,".json",sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    message("Loading cached file")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  config <- list()
  message(gatorURL)
  url <- gsub("/data/latest/", "/data/history/", gatorURL)
  if (!is.null(etag)) {
    config$httpheader['If-None-Match'] <- etag
  }
  message("Connecting to server to check for new data for ",fileId)

  file_request <- httr::GET(url,config=config)
  if (file_request$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return (origData)
  }
  if (file_request$status_code > 400 && file_request$status_code < 500) {
    message("Could not retrieved data from: ",url," got status code ",file_request$status_code)
    return ()
  }

  message(file_request$status_code)
  message("Retrieving data from Gator for ",httr::content(file_request)$title)
  retval <- rjson::fromJSON(httr::content(file_request,"text"))
  message("Retrieved data out from downloaded file")
  if ("etag" %in% names(retval)) {
    message("Setting etag retrieved from within response")
    retval$etag <- format(retval$etag,scientific=FALSE)
  } else {
    message("Setting etag retrieved from HTTP headers")
    retval$etag <- file_request$header[['etag']]
  }
  if (! "title" %in% names(retval)) {
    retval$title <- fileId
  }
  fileConn<-file(filename)
  writeLines(rjson::toJSON(retval), fileConn)
  close(fileConn)
  return (retval)
}

#' @importFrom rjson fromJSON
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom rjson toJSON
#' @importFrom httr add_headers
getGoogleFile <- function(fileId) {
  filename <- file.path(gator.cache,paste("gdrive-",fileId,".json",sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }

	access_info <- doSignin()
	gdrive_sig <- httr::add_headers(Authorization = paste('Bearer', access_info$access_token))
  url <- paste('https://www.googleapis.com/drive/v2/files/',fileId,sep='')
  if (!is.null(etag)) {
    gdrive_sig$httpheader['If-None-Match'] <- etag
  }
	file_meta <- httr::GET(url,config=gdrive_sig)
  if (file_meta$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return (origData)
  }
  if (file_meta$status_code > 400 && file_meta$status_code < 500) {
    message("Could not retrieve data from: ",url," got status code ",file_meta$status_code)
    return ()
  }
  message("Retrieving data from Google for ",httr::content(file_meta)$title)
	file_data <- httr::GET(httr::content(file_meta)$downloadUrl,gdrive_sig)
  retval <- rjson::fromJSON(httr::content(file_data,"text"))
  retval$etag <- httr::content(file_meta)$etag
  retval$title <- httr::content(file_meta)$title

  fileConn<-file(filename)
  writeLines(rjson::toJSON(retval), fileConn)
  close(fileConn)
	return (retval)
}

# Refresh OAuth 2.0 access token.
#
#' @importFrom httr POST
#' @importFrom httr content
oauth2.0_refresh <- function(endpoint, app, refresh_token, type = NULL) {
  # Use refresh_token to get a new (temporary) access token
  req <- httr::POST(
    url = endpoint$access,
    multipart = FALSE,
    body = list(
      client_id = app$key,
      client_secret = app$secret,
      grant_type = "refresh_token",
      refresh_token = refresh_token
    )
  )
  content_out <- httr::content(req, type = type)
  content_out <- c(
    content_out,
    refresh_token
  )
  stopifnot(
    length(content_out$expires_in) == 1
  )
  content_out$use_by <- Sys.time() + content_out$expires_in
  content_out
}
