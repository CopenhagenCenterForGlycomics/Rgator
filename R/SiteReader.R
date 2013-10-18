# We need to do this to load up the libraries to start with
# library(devtools)
# load_all('~/dev/Rgator/')

# Download to a DB with the etag / retrieved date?
# Store the key somewhere
# Make a parsing function attribute in the data file

#install_github("httr",username="hirenj")
library(httr)
library(plyr)
library(keychain)
library(RJSONIO)

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
	prefs = getPreferences()
	lapply(names(prefs$user_datasets),downloadDataset)
}

jsonParser <- function(data,keys) {
  currkeys <- unique(sapply(keys,FUN=function(key) { strsplit(key,".",fixed=TRUE)[[1]][1] }))
  keys <- unlist(lapply(keys,FUN=function(key) { sub(".*\\.","",key) }))
  if (length(currkeys) == 1) {
    return (ldply(data[[ currkeys[1] ]],.fun=function(dat) { parsed <- jsonParser(dat,keys); return(parsed); }))
  } else {
    datlist <- sapply(list(1:(length(keys)-1)),FUN=function(keyidx) {
      key <- keys[keyidx]
      if (is.null(data[[key]]) || length(data[[key]]) == 0) {
        data[[key]] = NULL
      }
      return (data[[key]])
    })
    key <- keys[length(keys)]
    if (is.null(data[[key]]) || length(data[[key]]) == 0) {
      data[[key]] = c(NULL)
    }
    datlist <- c( datlist, list(  data[[ keys[length(keys)]  ]] ))
    print (datlist)
    retval <- do.call("cbind",datlist)
    return (retval)
  }
}

downloadDataset <- function(set) {
  data <- getGoogleFile(set)
  assign(paste("gator.raw.",gsub("[[:space:]]|-","_",data$title),sep=""),data, envir = .GlobalEnv)
  if (!is.null(data$defaults$rKeys)) {
    frame <- ldply(data$data,.fun=function(dat) { return (jsonParser(dat,data$defaults$rKeys)) });
    names(frame) <- c("uniprot", data$defaults$rKeys)
  } else {
    frame <- ldply(data$data,.fun=function(dat) { return (ldply(dat$sites)) })
  }
  # sites

  # peptides.sequence,peptides.sites
  # function(dat) {
  #  ldply(dat$peptides,.fun=function(pep) {
  #    if (length(pep$sites) == 0) {
  #      pep$sites = NULL
  #    }
  #    return ( cbind(pep$sequence,pep$sites) )
  #  })
  # }
  currnames <- names(frame)
  currnames[1] <- 'uniprot'
  names(frame)<- currnames

  assign(paste("gator.",gsub("[[:space:]]|-","_",data$title),sep=""),frame, envir = .GlobalEnv)
  frame
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

	# access_info <- doSignin()
	# gdrive_sig <- sign_oauth2.0(access_info$access_token)
 #  url <- paste('https://www.googleapis.com/drive/v2/files/',fileId,sep='')
 #  if (!is.null(etag)) {
 #    gdrive_sig$httpheader['If-None-Match'] <- etag
 #  }
	# file_meta <- GET(url,config=gdrive_sig)
 #  if (file_meta$status_code == 304) {
    message("File data has not changed for ",origData$title)
    return (origData)
  # }
  message(file_meta$status_code)
  message("Retrieving data from Google for",content(file_meta)$title)
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