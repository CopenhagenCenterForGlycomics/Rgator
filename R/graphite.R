get_graphite_db <- function(db) {
  if (db == 'reactome') {
    return( graphite::reactome )
  }
  if (db == 'biocarta') {
    return( graphite::biocarta )
  }
  if (db == 'kegg') {
    return( graphite::kegg )
  }
  if (db == 'nci') {
    return( graphite::nci )
  }
  if (db == 'spike') {
    return( graphite::spike )
  }
  if (db == 'humancyc') {
    return( graphite::humancyc )
  }
  if (db == 'panther') {
    return( graphite::panther )
  }

}

#' @export
#' @importFrom graphite reactome biocarta kegg nci spike humancyc panther nodes
findCommonPathways.graphite <- function(organism=9606,max_pathway_size=30,...) {
  getBiocLiteLib('graphite')
  #require(graphite)
  genesets <- list(...)
  path.list <- sapply(genesets,function(genelist) {
    list <- subset(convertEntrezIds(organism,as.character(genelist)),uniprot %in% getUniprotIds(organism))
    list$uprotKey <- paste('UniProt:',toupper(list$uniprot),sep='')
    list$entrezKey <- paste('EntrezGene:',(list$geneid),sep='')
    list
  },USE.NAMES=T,simplify=F)
  for (name in names(path.list)) {
    attributes(path.list[[name]])$listname <- name
  }
  path_results <- rbindlist( lapply( c('reactome', 'biocarta', 'kegg','nci', 'spike', 'humancyc','panther' ), function(db) {
    res <- ldply(get_graphite_db(db),function(pw) {
      my.path.list <- path.list
      keyname <- "uprotKey"
      if (grepl("Entrez", graphite::nodes(pw)[1])) {
        keyname <- 'entrezKey'
      }
      if (length(graphite::nodes(pw) <= max_pathway_size)) {
        Reduce(function(left,right) {
          if (nrow(left) < 1) {
            return (left)
          }
          left$wantedkey <- left[[ keyname ]]
          right$wantedkey <- right[[ keyname ]]
          left_nodes <- subset( left, wantedkey %in% graphite::nodes(pw))
          right_nodes <- subset( right, wantedkey %in% graphite::nodes(pw))
          if (is.null(attributes(left)$listname) && nrow(right_nodes) > 0) {
            names(right_nodes) <- c(sapply( c('geneid','uniprot','genename','uprotKey','entrezKey'), function(x) { paste(x,attributes(right)$listname,sep='.')}),'wantedkey')
          }
          left <- merge(left_nodes,right_nodes,by='wantedkey',allow.cartesian=T,suffixes=c(paste('.',attributes(left)$listname,sep=''),paste('.',attributes(right)$listname,sep='')))
          if (nrow(left) > 0) {
            left[[keyname]] <- left$wantedkey
          }
          left$wantedkey <- NULL
          left$db <- rep(db,nrow(left))
          attributes(left)$listname <- NULL
          left
        }, my.path.list )
      }
    })
    return (res)
  }  ))
  names(path_results)[1] <- 'pathway'
  subset( path_results,select=which(! grepl("(uprot|entrez)Key",names(path_results)) ) )
}