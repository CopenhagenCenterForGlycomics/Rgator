# We need to do this to load up the libraries to start with
# library(devtools)
# load_all('~/dev/Rgator/')

#install_github("httr",username="hirenj")
library(httr)
library(plyr)
library(data.table)
library(RCurl)
library(rjson)
library(ggplot2)
library(VennDiagram)
library(gridExtra)


options(stringsAsFactors = FALSE)

current_token <- NULL
connection_key <- NULL

doSignin <- function() {
  if (!is.null(current_token)) {
    return((list(access_token=current_token)))
  }
  if (! is.null(getOption("current_token"))) {
    return(list(access_token=getOption("current_token")))
  }
	google <- oauth_endpoint(NULL,"auth","token", base_url= "https://accounts.google.com/o/oauth2")
	gdrive <- oauth_app("google",getOption("GoogleClientId"),secret=getOption("GoogleClientSecret"))

	token <- oauth2.0_token(google, gdrive,scope="https://www.googleapis.com/auth/drive.readonly",use_oob = TRUE )

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

askForSignin <- function(WS) {
  websocket_write(rjson::toJSON(list(message="upgradeConnection",data=getOption("connection_key"))), WS)
}

acceptToken <- function(json) {
  message("Successfully connected to GlycoDomain viewer")
  options(current_token = json$data$authtoken)
  options(connection_key = json$data$connectionkey)
}

gatorConnector <- function() {
  library(websockets)
  server = create_server(port=8880)
  daemonize(server)

  view <- function(prot) {
    if (length(server$client_sockets) > 0) {
      websocket_write(rjson::toJSON(list(message="showProtein", data=unique(prot))), server$client_sockets[[1]])
    }
  }

  receiver <- function(DATA,WS,...) {
    json <- rjson::fromJSON(rawToChar(DATA))
    if (json$message == "token") {
      acceptToken(json)
    }
  }

  setCallback("receive",receiver,server)
  setCallback("established", askForSignin, server)

  stopConnector <- function() {
    websocket_close(server)
    options(current_token = NULL)
    options(connection_key = NULL)
    rm("stopConnector",envir=.GlobalEnv)
    rm("Viewp",envir=.GlobalEnv)
  }
  assign("Viewp",view,envir = .GlobalEnv)
  assign("stopConnector",stopConnector,envir=.GlobalEnv)
}

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


jsonParser <- function(data,keys) {
  currkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key))>0) { strsplit(key,".",fixed=TRUE)[[1]][1] } }))
  currkeys <- unlist(currkeys[!sapply(currkeys,is.null)])

  localkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key)) < 1) { return(key) } }))
  localkeys <- unlist(localkeys[!sapply(localkeys,is.null)])

  keys <- unlist(lapply(keys,FUN=function(key) { sub("[^\\.]+\\.","",key) }))
  if (length(currkeys) == 1 && currkeys[1] != "*" && (length(keys) > length(currkeys) || keys[1] != currkeys[1]) ) {
    kidframe <- (ldply(data[[ currkeys[1] ]],.fun=function(dat) { parsed <- jsonParser(dat,keys[!keys %in% localkeys]); return(parsed); }))
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
        return (ldply(retval))
      }
      ##
      ## { "orthologs" : {"TAX1" : ["a","b","c"], "TAX2" : ["a","b"] }}
      ##
      if (length(unique(lapply(data[[keys[1]]],length))) > 1) {
        return (ldply(data[[keys[1]]],.fun=function(list) { return(data.frame(V1=list)) }))
      }
      return (ldply(data[[keys[1]]]))
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

downloadOrthologies <- function() {
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/homologene/',list(type='gatorURL',title='orthology.homologene'))
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/orthology.treefam.9/',list(type='gatorURL',title='orthology.treefam'))
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/orthology.inparanoid.8_1/',list(type='gatorURL',title='orthology.inparanoid'))
  assign('gator.orthology',rbind(gator.homologene,`gator.orthology.treefam`,`gator.orthology.inparanoid`),envir=.GlobalEnv)
}

downloadDomains <- function(organism) {
  downloadDataset('http://glycodomain-data.glycocode.com/data/latest/fulldomains/',list(type='gatorURL',title='fulldomains'))
  downloadDataset(paste('http://glycodomain-data.glycocode.com/data/latest/domains.',organism,'/',sep=''),list(type='gatorURL',title=paste('domains.',organism,sep='')))
}


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
            accumulated_frame <- rbindlist(list(accumulated_frame,frame))
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

  if (!is.null(data$defaults$rKeys) && length(data$defaults$rKeys) > 0) {
    all_prots <- names(data$data)
    frame <- rbindlist(llply(all_prots,.fun=function(uprot) {
      frm <- jsonParser(data$data[[uprot]],data$defaults$rKeys )

      # We should get a data frame out from the jsonParser - attach the uniprot id as
      # another column into the data frame

      frm$uniprot <- rep(uprot,dim(frm)[1])
      return(frm)
    },.progress="text"))

    # We need to re-arrange the columns here so that the uniprot column
    # ends up as the first column for consistency

    wanted_cols <- names(frame)
    frame <- frame[,c('uniprot',wanted_cols[!wanted_cols == 'uniprot'])]
    setnames(frame, c('uniprot', data$defaults$rKeys, rep('NA',dim(frame)[2] - (length(data$defaults$rKeys)+1))))
  } else {
    # Assume that we're just pulling out sites from the data sets if we're not given a particular
    # key to iterate over
    frame <- ldply(data$data,.fun=function(dat) { return (ldply(dat$sites)) },.progress="text")
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
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(basepath,showWarnings=FALSE)
  filename <- file.path(basepath,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  frame <- data.frame()

  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    origData <- load(fileConn)
    close(fileConn)
  }
  return(frame)
}

writeParsedJson <- function(frame,title) {
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(basepath,showWarnings=FALSE)
  filename <- file.path(basepath,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  save(frame,file=filename)
}

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
  frame <- rbindlist(llply(all_prots,.fun=function(uprot) {
    frm <- jsonParser(origData$data[[uprot]],origData$defaults$rKeys )

    # We should get a data frame out from the jsonParser - attach the uniprot id as
    # another column into the data frame

    frm$uniprot <- rep(uprot,dim(frm)[1])
    return(frm)
  },.progress="text"))
  # We need to re-arrange the columns here so that the uniprot column
  # ends up as the first column for consistency

  wanted_cols <- names(frame)
  frame <- frame[,c('uniprot',wanted_cols[!wanted_cols == 'uniprot'])]
  setnames(frame, c('uniprot', origData$defaults$rKeys, rep('NA',dim(frame)[2] - (length(origData$defaults$rKeys)+1))))

  if (!is.null(origData$defaults$rNames)) {
    names(frame) <- c('uniprot',origData$defaults$rNames)
  }

  return (frame)
}

getUniprotIds <- function(taxonomy) {
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
        accumulated_frame <- rbindlist(list(accumulated_frame,frame))
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
  fastas <- POST("http://www.uniprot.org/batch/",body=list(format='fasta',file=fileUpload('upload',toupper(paste(unlist(wanted_accs),collapse="\n")))),multipart=TRUE)
  if (fastas$status_code != 200) {
    message("Could not retrieve ids")
    return(unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
  }
  contents <- httr::content(fastas,"text")
  seqs <- strsplit(sub("\n","\t", unlist(strsplit(contents,"\n>"))),"\t")
  seqs <- ldply(seqs,function(row) { c(  row[1] , gsub("\n","",row[2]) )  });
  seqs$V1 <- sub(">?[st][pr]\\|","",seqs$V1)
  seqs$V1 <- sub("\\|.*","",seqs$V1)
  if (length(names(seqs)) < 2) {
    message("Could not retrieve sequences")
    message(accs)
    return (unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
  }
  names(seqs) <- c('uniprot','sequence')
  seqs$uniprot <- tolower(seqs$uniprot)
  assign('gator.UniProtData', rbindlist( list(get('gator.UniProtData'), seqs) ), envir = .GlobalEnv)
  writeParsedJson(gator.UniProtData,'UniProtData')
  Sys.sleep(wait)
  return (unique(subset(gator.UniProtData, uniprot %in% tolower(accs) )))
}

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

getDomainSets <- function( inputsites, sitecol, domaindata, max_dom_proportion=0.81, stem_distance=100  ) {
  message("Retrieving Uniprot sequences")
  seqdat <- getUniprotSequences(unique(inputsites$uniprot))
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
  typeii <- unique(ddply(between,.(sitekey),function(input) {
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
  stem_typeii <- ddply(between,.(sitekey),function(input) {
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
  signalp_stem <- ddply(between,.(sitekey),function(input) {
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
  stem_typei <- ddply(subset(between, ! sitekey %in% signalp_stem$sitekey),.(sitekey),function(input) {
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

uniqueframe <- function(set){
  return(as.data.frame(unique(as.matrix(set))))
}

generateLogoPlot <- function(dataframe,windowcol) {
  uniprot_2013_12_freq <- list(A=0.0825,R=0.0553,N=0.0406,D=0.0545,C=0.0137,Q=0.0393,E=0.0675,G=0.0707,H=0.0227,I=0.0595,L=0.0966,K=0.0584,M=0.0242,F=0.0386,P=0.0470,S=0.0657,T=0.0534,W=0.0108,Y=0.0292,V=0.0686)
  pwm <- calculatePWM(dataframe,windowcol,names(uniprot_2013_12_freq))
  return(berrylogo(pwm,uniprot_2013_12_freq))
}

generateVennDiagram <- function(data=list(),title="Venn Diagram") {
  return (gTree(children = gList(grid.grabExpr(grid.draw(venn.diagram(data,filename=NULL,main=title)))), cl=c("arrange", "ggplot")))
}

berrylogo<-function(pwm,backFreq,zero=.0001){
  pwm[pwm==0]<-zero
  bval<-plyr::laply(names(backFreq),function(x){log(pwm[x,])-log(backFreq[[x]])})
  row.names(bval)<-names(backFreq)
  window_size = floor( 0.5*length(dimnames(bval)[[2]]) )
  dimnames(bval)[[2]]<- c((-1*window_size):window_size)
  p<-ggplot2::ggplot(reshape2::melt(bval,varnames=c("aa","pos")),ggplot2::aes(x=pos,y=value,label=aa))+
    ggplot2::geom_abline(ggplot2::aes(slope=0), colour = "grey",size=2)+
    ggplot2::geom_text(ggplot2::aes(colour=factor(aa)),size=8)+
    ggplot2::theme(legend.position="none")+
    ggplot2::scale_x_continuous(name="Position",breaks=(-1*window_size):window_size)+
    ggplot2::scale_y_continuous(name="Log relative frequency")
  return(p)
}

getEntrezIds <- function(organism,ids) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dbname<-organisms[[as.character(organism)]]
  if ( ! library(dbname,lib.loc=c(basepath),character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    biocLite(dbname,lib=basepath)
  }
  library(dbname,character.only=TRUE,lib.loc=c(basepath))
  uprotmap <- sub("\\.db","UNIPROT",dbname)
  uprots <- unique(intersect(  toupper(ids),mappedRkeys(get(uprotmap))))
  entrez_ids <- toTable(revmap(get(uprotmap))[  uprots ])
  names(entrez_ids) <- c('gene_id','uniprot')
  gene_ids <- entrez_ids$gene_id
  names(gene_ids) <- entrez_ids$uniprot
  return (gene_ids)
}

convertEntrezIds <- function(organism,ids) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dbname<-organisms[[as.character(organism)]]
  if ( ! library(dbname,lib.loc=c(basepath),character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    biocLite(dbname,lib=basepath)
  }
  library(dbname,character.only=TRUE,lib.loc=c(basepath))
  retdata <- select(get(dbname), keys=as.character(ids), cols=c('UNIPROT', 'SYMBOL', 'ENTREZID'), keytype="ENTREZID")
  names(retdata) <- c('geneid','uniprot','genename')
  retdata
}

getGOenrichmentGenes <- function(enrichment,wanted_terms=c(),organism=9606) {
  unique(convertEntrezIds(organism,unique(unlist(sapply(wanted_terms, function(go) {   geneIdsByCategory(enrichment)[[go]]  })))))
}

getGOTerms <- function(organism,uniprots) {
  if ( ! library("GO.db",character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    biocLite("GO.db")
  }
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  query_ids <- getEntrezIds(organism,uniprots)
  godb <- get( sub("\\.db","GO",organisms[as.character(organism)] ) )
  terms <- toTable(godb[query_ids])
  names(terms) <- c('gene_id','go_id','Evidence','Ontology')
  terms <- merge( terms, data.frame(gene_id=query_ids,uniprot=tolower(names(query_ids))), by='gene_id')
  library('GO.db')
  if (dim(terms)[1] > 0) {
    go_terms <- (select(GO.db, terms$go_id,"TERM"))
    names(go_terms) <- c('go_id','term')
    terms <- merge( terms, go_terms, by='go_id')
    return (terms)
  }
  return (terms)
}

GO_children <- function(node = "GO:0008150", ontology = "BP") {
    library('GO.db')
    if (ontology == "BP") GOCHILDREN <- GOBPCHILDREN
    if (ontology == "CC") GOCHILDREN <- GOCCCHILDREN
    if (ontology == "MF") GOCHILDREN <- GOMFCHILDREN
    parents <- node

    # initialize output
    out <- c(parents)

    # do the following until there are no more parents
    while (any(!is.na(parents))) {
        # Get the unique children of the parents (that aren't NA)
        children <- unique(unlist(mget(parents[!is.na(parents)], GOCHILDREN)))

        # append chldren to beginning of `out`
        # unique will keep the first instance of a duplicate
        # (i.e. the most recent child is kept)
        out <- unique(append(children[!is.na(children)], out))

        # children become the parents of the next generation
        parents <- children
    }
    return(out)
}

getGOEnrichment <- function(organism,uniprots,query_ids=c(),universe=c(),ontology='BP',direction='over') {
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  if ( ! library("GO.db",character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    biocLite("GO.db")
  }
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  if ( ! library(dbname,lib.loc=c(basepath),character.only=TRUE,logical.return=TRUE,quietly=TRUE)) {
    biocLite(dbname,lib=basepath)
  }
  library(dbname,character.only=TRUE,lib.loc=c(basepath))
  if (length(query_ids) < 1 & length(uniprots) > 0) {
    query_ids <- getEntrezIds(organism,uniprots)
  }
  if (length(universe) < 1) {
    universe <- as.list( mappedkeys(get( sub("\\.db","GO",organisms[as.character(organism)] ) )) )
  } else {
    universe <- getEntrezIds(organism,universe)
  }
  library('GO.db')
  library('GOstats')
  params <- new('GOHyperGParams',
              geneIds=query_ids,
              universeGeneIds=unlist(universe),
              ontology=ontology,
              pvalueCutoff=0.05,
              conditional=TRUE,
              testDirection=direction,
              annotation=as.character(organisms[as.character(organism)])
             )
  hgOver <- hyperGTest(params)
  return (hgOver)
}


addSiteColumn <- function(dataframe) {
  dataset <- ddply(dataframe,.(uniprot),function(df) { vals <- c(1:(dim(df)[1])); df$site <- vals; return (df); })
  eval.parent(substitute(dataframe<-dataset))
}

getGatorSnapshotSubset <- function(fileId,accs) {
  url <- paste('http://localhost:3001/data/history/',fileId,'?accs=',paste(tolower(unlist(accs)),collapse=','),sep='')
  config <- list()
  file_request <- GET(url,config=config)
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
  retval <- httr::content(file_request)
  retval$etag <- format(retval$etag,scientific=FALSE)
  return (retval)
}

getGatorSnapshot <- function(gatorURL,fileId) {
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(basepath,showWarnings=FALSE)
  filename <- file.path(basepath,paste("gator-",fileId,".json",sep=''))
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

  file_request <- GET(url,config=config)
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
  retval <- httr::content(file_request)
  message("Retrieved data out from downloaded file")
  if ("etag" %in% names(retval)) {
    message("Setting etag retrieved from within response")
    retval$etag <- format(retval$etag,scientific=FALSE)
  } else {
    message("Setting etag retrieved from HTTP headers")
    retval$etag <- file_request$header[['etag']]
  }
  if (! "title" %in% names(retval) && title ) {
    retval$title <- fileId
  }
  fileConn<-file(filename)
  writeLines(rjson::toJSON(retval), fileConn)
  close(fileConn)
  return (retval)
}


getGoogleFile <- function(fileId) {
  basepath <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(basepath,showWarnings=FALSE)
  filename <- file.path(basepath,paste("gdrive-",fileId,".json",sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }

	access_info <- doSignin()
	gdrive_sig <- add_headers(Authorization = paste('Bearer', access_info$access_token))
  url <- paste('https://www.googleapis.com/drive/v2/files/',fileId,sep='')
  if (!is.null(etag)) {
    gdrive_sig$httpheader['If-None-Match'] <- etag
  }
	file_meta <- GET(url,config=gdrive_sig)
  if (file_meta$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return (origData)
  }
  if (file_meta$status_code > 400 && file_meta$status_code < 500) {
    message("Could not retrieve data from: ",url," got status code ",file_meta$status_code)
    return ()
  }
  message("Retrieving data from Google for ",httr::content(file_meta)$title)
	file_data <- GET(httr::content(file_meta)$downloadUrl,gdrive_sig)
  retval <- httr::content(file_data)
  retval$etag <- httr::content(file_meta)$etag
  retval$title <- httr::content(file_meta)$title

  fileConn<-file(filename)
  writeLines(toJSON(retval), fileConn)
  close(fileConn)
	return (retval)
}

# Refresh OAuth 2.0 access token.
#
oauth2.0_refresh <- function(endpoint, app, refresh_token, type = NULL) {
  # Use refresh_token to get a new (temporary) access token
  req <- POST(
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
