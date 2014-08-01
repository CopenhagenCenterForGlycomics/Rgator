
graphite.reactome <- function(organism=9606,max_pathway_size=30,...) {
  getBiocLiteLib('graphite')
  #require(graphite)
  genesets <- list(...)
  path.list <- sapply(genesets,function(genelist) {
    list <- subset(convertEntrezIds(organism,as.character(genelist)),uniprot %in% getUniprotIds(organism))
    list$uprotKey <- paste('UniProt:',toupper(list$uniprot),sep='')
    list$entrezKey <- paste('EntrezGene:',(list$geneid),sep='')
    list
  },USE.NAMES=T,simplify=F)
  path_results <- rbindlist( lapply( c('reactome', 'biocarta', 'kegg','nci', 'spike', 'humancyc','panther' ), function(db) {
    res <- ldply(get(db,getNamespace('graphite')),function(pw) {
      my.path.list <- path.list
      keyname <- "uprotKey"
      if (grepl("Entrez", nodes(pw)[1])) {
        keyname <- 'entrezKey'
      }
      if (length(nodes(pw) <= max_pathway_size)) {
        Reduce(function(left,right) {
          my.path.list[[left]]$wantedkey <- my.path.list[[left]][[ keyname ]]
          my.path.list[[right]]$wantedkey <- my.path.list[[right]][[ keyname ]]

          left_nodes <- subset( my.path.list[[left]], wantedkey %in% nodes(pw))
          right_nodes <- subset( my.path.list[[right]], wantedkey %in% nodes(pw))
          my.path.list[[left]] <- merge(left_nodes,right_nodes,by='wantedkey',allow.cartesian=T,suffixes=c(paste('.',left,sep=''),paste('.',right,sep='')))
          my.path.list[[left]]$wantedkey <- NULL
          my.path.list[[left]]$db <- rep(db,nrow(my.path.list[[left]]))
          my.path.list[[left]] <- my.path.list[[left]]
        }, names(my.path.list) )
      }
    })
    return (res)
  }  ))
  names(path_results)[1] <- 'pathway'
  subset( path_results,select=which(! grepl("Key\\.",names(path_results))) )
}