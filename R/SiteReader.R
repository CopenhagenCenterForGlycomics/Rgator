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

# @importFrom httr oauth_endpoint
# @importFrom httr oauth_app
# @importFrom httr oauth2.0_token
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

# @importFrom rjson fromJSON
# @importFrom websockets websocket_write
askForSignin <- function(socket) {
  socket$send(rjson::toJSON(list(message="upgradeConnection",data=getOption("connection_key"))))
}

acceptToken <- function(json) {
  message("Successfully connected to GlycoDomain viewer")
  options(current_token = json$data$authtoken)
  options(connection_key = json$data$connectionkey)
}

acceptPreferences <- function(json) {
  if (json$data$connectionkey == getOption("connection_key")) {
    fileConn<-file('auto.domaintoolsession')
    writeLines(rjson::toJSON(json$data$preferences), fileConn)
    close(fileConn)
    message("Successfully received preferences from GlycoDomain viewer")
  }
}

# @importFrom websockets daemonize
# @importFrom websockets create_server
# @importFrom websockets websocket_write
# @importFrom websockets websocket_close
# @importFrom rjson fromJSON
# @importFrom websockets setCallback
#' @export
gatorConnector <- function() {
  assign('gator.sockets',list(), envir = .GlobalEnv)
  server <- httpuv::startDaemonizedServer(host='0.0.0.0',port=8880,list(onWSOpen=function(ws) {
    gator.sockets[[length(gator.sockets) + 1]] <- ws
    socket_id <- length(gator.sockets)
    assign('gator.sockets',gator.sockets, envir = .GlobalEnv)
    ws$onMessage(function(binary,data) {
      receiver(data)
    })
    askForSignin(ws)
    ws$onClose(function() {
      gator.sockets[[socket_id]] <- NULL
      assign('gator.sockets',gator.sockets, envir = .GlobalEnv)
    })
  }))

  view <- function(prot) {
    message(length(gator.sockets))
    for (sock in gator.sockets) {
      if (! is.null(sock)) {
        sock$send(rjson::toJSON(list(message="showProtein", data=unique(prot))))
      }
    }
  }
  align <- function(prot) {
    for (sock in gator.sockets) {
      if (! is.null(sock)) {
        sock$send(rjson::toJSON(list(message="alignProtein", data=unique(prot))))
      }
    }
  }
  compactRenderer <- function() {
    for (sock in gator.sockets) {
      if (! is.null(sock)) {
        sock$send(rjson::toJSON(list(message="compactRenderer", data="")))
      }
    }
  }

  getPreferences <- function() {
    if (length(gator.sockets) > 0) {
        sock <- gator.sockets[[1]]
        if (! is.null(sock)) {
          sock$send(rjson::toJSON(list(message="retrieveSession", data=getOption("connection_key"))))
        }
    }
  }

  receiver <- function(DATA) {
    json <- rjson::fromJSON(DATA)
    if (json$message == "token") {
      acceptToken(json)
    }
    if (json$message == "preferences") {
      acceptPreferences(json)
    }
  }


  stopConnector <- function() {
    httpuv::stopDaemonizedServer(server)
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
  assign("getPreferences",getPreferences,envir=.GlobalEnv)
}


# @importFrom rjson fromJSON
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

  updateDataVersions()

}

# @importFrom plyr ldply
jsonParser <- function(data,keys) {
  options(stringsAsFactors = FALSE)
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
  data.env = getDataEnvironment()
  data.env[[ 'gator.orthology' ]] <- rbind(gator.homologene,`gator.orthology.treefam`,`gator.orthology.inparanoid`)
}

has_internet <- function() {
  url = 'http://www.google.com'
  test <- try(suppressWarnings(readLines(url,n=1)),silent=TRUE)

  !inherits(test,'try-error')
}

# @importFrom plyr ldply
# @importFrom plyr llply
# @importFrom data.table rbindlist
# @importFrom data.table setnames
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
            data.env = getDataEnvironment()
            data.env[[ config[['title']] ]] <- accumulated_frame
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

  if (is.null(data)) {
    data<-list(etag=NULL,title=gsub("[[:space:]]|-","_",config[['title']]))
  }

  # We should check to make sure that the etag for the current data frame in the global namespace
  # matches with the etag for the current data, so we can find out if we have to redo the parsing of
  # the data
  message("Checking etags of pre-parsed data")
  if(etagcheck) {
    if (exists( data$title, where=getDataEnvironment() )) {
      message("We have data loaded up in the environment")
      frame <- get( data$title, envir=getDataEnvironment() )
      if ( is.null(data$etag) ) {
        return ()
      }
      etagmatch=!is.null(attributes(frame)$etag) && attributes(frame)$etag == format(data$etag,scientific=FALSE)
      if ( is.null(data$etag) || etagmatch ) {
        return ()
      }
      message("No etag match, continuing")
    } else {
      loadParsedJson(data$title)
      frame <- getDataEnvironment()[[data$title]]
      etagmatch=!is.null(attributes(frame)$etag) && attributes(frame)$etag == format(data$etag,scientific=FALSE)
      if ( is.null(data$etag) || etagmatch ) {
        message("We have data that has already been parsed")
        data.env = getDataEnvironment()
        data.env[[ data$title ]] <- frame
        return ()
      }

    }
  }

  message("Preparing to parse data")

  parserFunction = jsonParser
  version = NULL

  if ('metadata' %in% names(data) && 'msdata-version' %in% data$metadata) {
    message("Parsing msdata")
    parserFunction = chooseMsDataParser(data$metadata[['msdata-version']])
    version = getMsDataVersionId(data$metadata[['msdata-version']],data)
  }
  if ('metadata' %in% names(data) && is.list(data$metadata) && 'msdata-version' %in% names(data$metadata[[1]])) {
    message("Parsing msdata")
    parserFunction = chooseMsDataParser(data$metadata[[1]][['msdata-version']])
    version = getMsDataVersionId(data$metadata[[1]][['msdata-version']],data)
  }

  if (!is.null(data$defaults$rKeys) && length(data$defaults$rKeys) == 1 && 'base64' %in% data$defaults$rKeys) {
    all_prots <- names(data$data)
    names(all_prots) <- all_prots
    frame <- plyr::llply(all_prots,.fun=function(uprot) {
      if (head(nchar(data$data[[uprot]]),1) > 0) {
        decodeBase64(data$data[[uprot]])
      }
    },.progress="text")
    attributes(frame)$etag <- data$etag
    data.env = getDataEnvironment()
    data.env[[ filename ]] <- frame
    if ( ! is.null(version)) {
      atttributes(frame)$version <- version
    }
    writeParsedJson(filename)
    return (data.frame())
  } else if ( (!is.null(data$defaults$rKeys) && length(data$defaults$rKeys) > 0) || ('metadata' %in% names(data)) ) {
    all_prots <- names(data$data)
    frame <- data.table::rbindlist(plyr::llply(all_prots,.fun=function(uprot) {
      frm <- parserFunction(data$data[[uprot]],data$defaults$rKeys )
      # We should get a data frame out from the jsonParser - attach the uniprot id as
      # another column into the data frame

      frm$uniprot <- rep(uprot,dim(frm)[1])
      if (nrow(frm) < 1) {
        return (NULL)
      }
      return(frm)
    },.progress="text"))

    # We need to re-arrange the columns here so that the uniprot column
    # ends up as the first column for consistency
    wanted_cols <- names(frame)
    frame <- subset(frame,select=c('uniprot',wanted_cols[!wanted_cols == 'uniprot']))
    if ("defaults" %in% names(data) && ! is.null(data$defaults$rKeys)) {
      data.table::setnames(frame, c('uniprot', data$defaults$rKeys, rep('NA',dim(frame)[2] - (length(data$defaults$rKeys)+1))))
    }
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

  frame$uniprot <- tolower(frame$uniprot)

  if ( ! is.null(version) ) {
    attributes(frame)$version <- version
  }

  data.env = getDataEnvironment()
  data.env[[ data$title ]] <- frame

  writeParsedJson(data$title)

  frame
}

loadParsedJson <- function(title) {
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  frame <- data.frame(stringsAsFactors=F)

  if (file.exists(filename)) {
    origData <- readRDS(filename)
    data.env = getDataEnvironment()
    data.env[[ title ]] <- origData
  }
  return(frame)
}

writeParsedJson <- function(title) {
  message("Writing out cache for ",title)
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  saveRDS(get(title,envir=as.environment('package:gatordata')),filename)
}

# @importFrom RCurl base64Decode
decodeBase64 <- function(base64) {
  decoded <- RCurl::base64Decode(base64,'raw');
  readBin( decoded, 'numeric', length(decoded)/4 ,4);
}

# @importFrom rjson fromJSON
# @importFrom plyr llply
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
  frame <- plyr::llply(all_prots,.fun=function(uprot) {
    if (head(nchar(origData$data[[uprot]]),1) > 0) {
      if (length(origData$data[[uprot]]) > 1) {
        plyr::llply( origData$data[[uprot]], decodeBase64 )
      } else {
        decodeBase64(origData$data[[uprot]])
      }
    }
  },.progress="text")

  return (frame)
}

# @importFrom rjson fromJSON
# @importFrom data.table rbindlist
# @importFrom data.table setnames
# @importFrom plyr llply
testParseJson <- function(filename,attach=F) {
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    message("Loading cached file")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  all_prots <- names(origData$data)
  parserFunction = jsonParser

  if ( 'metadata' %in% names(origData) && 'msdata-version' %in% names(origData$metadata[[1]]) ) {
    parserFunction = chooseMsDataParser(origData$metadata[[1]][['msdata-version']])
  }

  frame <- data.table::rbindlist(plyr::llply(all_prots,.fun=function(uprot) {
    frm <- parserFunction(origData$data[[uprot]],origData$defaults$rKeys )

    # We should get a data frame out from the jsonParser - attach the uniprot id as
    # another column into the data frame

    frm$uniprot <- rep(uprot,dim(frm)[1])
    return(frm)
  },.progress="text"))

  # We need to re-arrange the columns here so that the uniprot column
  # ends up as the first column for consistency

  wanted_cols <- names(frame)
  frame <- subset(frame,select=c('uniprot',wanted_cols[!wanted_cols == 'uniprot']))
  if( 'defaults' %in% names(origData) && 'rKeys' %in% names(origData$defaults) ) {
    data.table::setnames(frame, c('uniprot', origData$defaults$rKeys, rep('NA',dim(frame)[2] - (length(origData$defaults$rKeys)+1))))
  }

  if (!is.null(origData$defaults$rNames)) {
    names(frame) <- c('uniprot',origData$defaults$rNames)
  }
  frame <- as.data.frame(frame)

  if (attach) {
    data.env = getDataEnvironment()
    data.env[[ filename ]] <- frame
  }

  return (as.data.frame(frame))
}

getDataEnvironment <- function() {
  if (! 'package:gatordata' %in% search()) {
    attach(new.env(),name='package:gatordata')
    dataenv <- as.environment('package:gatordata')
    package_path <- file.path(tempdir(), 'gatordata-package')
    dir.create(package_path, showWarnings = FALSE)
    attr(dataenv,'path') <- package_path
    updateDataVersions()
  }
  return(as.environment('package:gatordata'))
}

getDataVersions <- function() {
  env <- getDataEnvironment()
  datasets <- ls(env)
  version_ids <- sapply(datasets,function(set) {
    if ('version' %in% names(attributes(env[[set]]))) {
      paste( set, attributes(env[[set]])[['version']], sep='-' )
    } else {
      NA
    }
  },USE.NAMES=F)
  version_ids[ ! is.na(version_ids) ]
}

updateDataVersions <- function() {
  package_path <- attr(getDataEnvironment(),'path')
  description <- list(Package='gatordata',Type='Package',Title='gatordata',Version=paste(getDataVersions(),collapse='_'),Date=format(Sys.time(), "%y-%m-%d"))
  write(renderDescription(description),file.path(package_path,'DESCRIPTION'))
}

renderDescription <- function(fields) {
  paste(sapply(names(fields),function(x){ paste(x,fields[[x]],sep=':'); },USE.NAMES=F),collapse="\n")
}

cacheFile <- function(url,fileId,gzip=F,...) {
  filename <- file.path(gator.cache,paste("gator-",fileId,sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    if (gzip) {
      filename <- gzfile(filename)
    }
    return (read.delim(filename,header=F,sep='\t',...))
  }
  download.file(url,filename)

  if (gzip) {
    filename <- gzfile(filename)
  }
  return (read.delim(filename,header=F,sep='\t',...))
}

uniqueframe <- function(set){
  return(as.data.frame(unique(as.matrix(set))))
}

getBiocLiteLib <- function(dbname) {
  if ( ! suppressWarnings(require(dbname,character.only=TRUE,quietly=TRUE))) {
    if ( gator.biocLite.stdlib) {
      biocLite(dbname)
      return
    }
  }
  if ( ! suppressWarnings(require(dbname,lib.loc=c(gator.cache),character.only=TRUE,quietly=TRUE)))  {
    if (! gator.biocLite.stdlib) {
      biocLite(dbname,lib=gator.cache)
      return
    }
  }

}



# @importFrom plyr ddply
# @importFrom plyr .
addSiteColumn <- function(dataframe) {
  dataset <- plyr::ddply(dataframe,plyr::.(uniprot),function(df) { vals <- c(1:(nrow(df))); df$site <- vals; return (df); })
  eval.parent(substitute(dataframe<-dataset))
}

# @importFrom rjson fromJSON
# @importFrom httr GET
# @importFrom httr content
getGatorSnapshotSubset <- function(fileId,accs) {
  url <- paste('http://localhost:3001/data/history/',fileId,'?accs=',paste(tolower(unlist(accs)),collapse=','),sep='')
  config <- list()
  file_request <- httr::GET(url,config=config)
  if (file_request$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return ()
  }
  if (file_request$status_code > 400 && file_request$status_code < 500) {
    message("Could not retrieve data from: ",url," got status code ",file_request$status_code)
    return ()
  }

  message(file_request$status_code)
  message("Retrieving data from Gator for ",httr::content(file_request)$title)
  retval <- rjson::fromJSON(httr::content(file_request,"text"))
  retval$etag <- format(retval$etag,scientific=FALSE)
  return (retval)
}

# @importFrom rjson fromJSON
# @importFrom httr GET
# @importFrom httr content
getGatorSnapshot <- function(gatorURL,fileId) {
  filename <- file.path(gator.cache,paste("gator-",fileId,".json",sep=''))
  etag <- NULL
  origData <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    message("Loading cached file")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  if (! has_internet()) {
    return (origData)
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
    message("Could not retrieve data from: ",url," got status code ",file_request$status_code)
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

# @importFrom rjson fromJSON
# @importFrom httr GET
# @importFrom httr content
# @importFrom rjson toJSON
# @importFrom httr add_headers
getGoogleFile <- function(fileId) {
  filename <- file.path(gator.cache,paste("gdrive-",fileId,".json",sep=''))
  etag <- NULL
  if (file.exists(filename)) {
    fileConn <- file(filename,"r")
    origData <- rjson::fromJSON(,fileConn)
    close(fileConn)
    etag <- origData$etag
  }
  if (! has_internet()) {
    return (origData)
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
# @importFrom httr POST
# @importFrom httr content
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
