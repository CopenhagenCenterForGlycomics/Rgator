# We need to do this to load up the libraries to start with
# library(devtools)
# load_all('~/dev/Rgator/')

getBiocLiteLib <- function(dbname) {
  if ( ! suppressWarnings(require(dbname,character.only=TRUE,quietly=TRUE))) {
    biocLite(dbname)
  }
}
