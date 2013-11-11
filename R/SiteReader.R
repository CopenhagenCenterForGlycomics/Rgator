# We need to do this to load up the libraries to start with
# library(devtools)
# load_all('~/dev/Rgator/')

#install_github("httr",username="hirenj")
library(httr)
library(plyr)
library(keychain)
library(RJSONIO)
library(data.table)
library(RCurl)

options(stringsAsFactors = FALSE)

current_token <- NULL

doSignin <- function() {
  if (!is.null(current_token)) {
    return(current_token)
  }
	google <- oauth_endpoint(NULL,"auth","token", base_url= "https://accounts.google.com/o/oauth2")
	gdrive <- oauth_app("google",getOption("GoogleClientId"),secret=getOption("GoogleClientSecret"))

	if (!is.null(getPassword("GoogleRefreshToken",quiet=TRUE))) {
		return (oauth2.0_refresh(google,gdrive,getPassword("GoogleRefreshToken")))
	}
	token <- oauth2.0_token(google, gdrive,scope="https://www.googleapis.com/auth/drive.readonly",use_oob = TRUE )
	setPassword("GoogleRefreshToken",,token$refresh_token)
  current_token <- token
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

syncDatasets <- function() {
	prefs = getPreferences('0By48KKDu9leCVEFVdS0xOVNSblE')
	lapply(seq_along(prefs$user_datasets), function(i) { key <- names(prefs$user_datasets)[i]; config <- prefs$user_datasets[[key]]; downloadDataset(key,config) } )
}

jsonParser <- function(data,keys) {
  currkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key))>0) { strsplit(key,".",fixed=TRUE)[[1]][1] } }))
  currkeys <- unlist(currkeys[!sapply(currkeys,is.null)])

  localkeys <- unique(sapply(keys,FUN=function(key) { if (length(grep("\\.",key)) < 1) { return(key) } }))
  localkeys <- unlist(localkeys[!sapply(localkeys,is.null)])

  keys <- unlist(lapply(keys,FUN=function(key) { sub("[^\\.]+\\.","",key) }))
  if (length(currkeys) == 1 && currkeys[1] != "*" && (length(keys) > length(currkeys) || keys[1] != currkeys[1]) ) {
    kidframe <- (ldply(data[[ currkeys[1] ]],.fun=function(dat) { parsed <- jsonParser(dat,keys[!keys %in% localkeys]); return(parsed); }))
    for (local in localkeys) {
      kidframe[[local]] <- rep(data[[local]],dim(kidframe)[1])
    }
    return (kidframe)
  } else {
    if (length(keys) == 1) {
      if (length(currkeys) > 0 && currkeys[1] == "*") {
        retval <-lapply( data, FUN=function(dat) { return ((jsonParser(dat,keys))) }  )
        return (ldply(retval))
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

downloadDataset <- function(set,config,accs=c(),etagcheck=TRUE) {
  if (length(grep("gatorURL",config[['type']])) > 0) {
    # This is only the key for the hash in the prefs file
    # We need to be storing the type with the key
    # set.type == gatorURL
    # If there's a trailing slash - we want to do a query
    if (substr(set, nchar(set), nchar(set)) == '/' || length(accs) > 0) {
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
      data <- getGatorSnapshot(set,gsub("[[:space:]]|-","_",config[['title']]))
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

  assign(paste("gator.raw.",gsub("[[:space:]]|-","_",data$title),sep=""),data, envir = .GlobalEnv)

  # We should check to make sure that the etag for the current data frame in the global namespace
  # matches with the etag for the current data, so we can find out if we have to redo the parsing of
  # the data

  if(etagcheck && exists( paste("gator.",gsub("[[:space:]]|-","_",data$title),sep="") )) {
    frame <- get( paste("gator.",gsub("[[:space:]]|-","_",data$title),sep="") )
    if (!is.null(attributes(frame)$etag) && attributes(frame)$etag == format(data$etag,scientific=FALSE)) {
      return ()
    }
  }


  if (!is.null(data$defaults$rKeys) && length(data$defaults$rKeys) > 0) {
    all_prots <- names(data$data)
    frame <- ldply(all_prots,.fun=function(uprot) {
      frm <- jsonParser(data$data[[uprot]],data$defaults$rKeys )

      # We should get a data frame out from the jsonParser - attach the uniprot id as
      # another column into the data frame

      frm$uniprot <- rep(uprot,dim(frm)[1])
      return(frm)
    },.progress="text")

    # We need to re-arrange the columns here so that the uniprot column
    # ends up as the first column for consistency

    wanted_cols <- names(frame)
    frame <- frame[,c('uniprot',wanted_cols[!wanted_cols == 'uniprot'])]
    names(frame) <- c('uniprot', data$defaults$rKeys, rep(NA,dim(frame)[2] - (length(data$defaults$rKeys)+1)))
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
  frame
}

getUniprotSequences <- function(accs) {
  wanted_accs <- accs
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
  contents <- content(fastas)
  seqs <- strsplit(sub("\n","\t", unlist(strsplit(contents,"\n>"))),"\t")
  seqs <- ldply(seqs,function(row) { c(  row[1] , gsub("\n","",row[2]) )  });
  seqs$V1 <- sub(">?sp\\|","",seqs$V1)
  seqs$V1 <- sub("\\|.*","",seqs$V1)
  names(seqs) <- c('uniprot','sequence')
  assign('gator.UniProtData', rbindlist( list(get('gator.UniProtData'), seqs) ), envir = .GlobalEnv)
  return (subset(gator.UniProtData, uniprot %in% accs ))
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
  message("Retrieving data from Gator for ",content(file_request)$title)
  retval <- content(file_request)
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
    origData <- fromJSON(fileConn,simplifyWithNames=FALSE)
    close(fileConn)
    etag <- origData$etag
  }
  config <- list()
  message(gatorURL)
  url <- gsub("/data/latest/", "/data/history/", gatorURL)
  if (!is.null(etag)) {
    config$httpheader['If-None-Match'] <- etag
  }
  message("Connecting to server to retrieve data for ",fileId)

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
  message("Retrieving data from Gator for ",content(file_request)$title)
  retval <- content(file_request)
  message("Retrieved data out from downloaded file")
  retval$etag <- format(retval$etag,scientific=FALSE)
  fileConn<-file(filename)
  writeLines(toJSON(retval), fileConn)
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
    origData <- fromJSON(fileConn,simplifyWithNames=FALSE)
    close(fileConn)
    etag <- origData$etag
  }

	access_info <- doSignin()
	gdrive_sig <- sign_oauth2.0(access_info$access_token)
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
    message("Could not retrieved data from: ",url," got status code ",file_meta$status_code)
    return ()
  }
  message("Retrieving data from Google for ",content(file_meta)$title)
	file_data <- GET(content(file_meta)$downloadUrl,gdrive_sig)
  retval <- content(file_data)
  retval$etag <- content(file_meta)$etag
  retval$title <- content(file_meta)$title

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
  content_out <- content(req, type = type)
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