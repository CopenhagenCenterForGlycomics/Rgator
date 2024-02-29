loadParsedJson <- function(title) {
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  message("Reading cache for ",title," from ",filename)
  frame <- data.frame(stringsAsFactors=F)

  if (file.exists(filename)) {
    origData <- readRDS(filename)
    data.env = getDataEnvironment()
    data.env[[ title ]] <- origData
  }
  return(frame)
}

writeParsedJson <- function(title) {
  filename <- file.path(gator.cache,paste("gator.parsed.",gsub("[[:space:]]|-","_",title),sep=""))
  message("Writing out cache for ",title," to ",filename)
  saveRDS(get(title,envir=getDataEnvironment()),filename)
}

loaded_data = new.env()

getDataEnvironment <- function() {
  loaded_data
}

gator.cache <- getOption("gator.cache")

if (! is.null(gator.cache)) {
  dir.create(gator.cache,showWarnings=FALSE)
  if ( ! gator.cache %in% .libPaths() ) {
    .libPaths(c(.libPaths(), gator.cache ))
  }
} else {
  gator.cache <- file.path(system.file(package="Rgator"),"cachedData")
  dir.create(gator.cache,showWarnings=FALSE)
}