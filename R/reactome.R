#' @export
findCommonPathways <- function(organism=9606,max_pathway_size=80,...) {
  rbind( findCommonPathways.reactome(organism,max_pathway_size,...), findCommonPathways.graphite(organism,max_pathway_size,...) )
}

reactome.total_count <- function(pathid) {
  nrow(AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID[pathid]))
}

#' @export
findCommonPathways.reactome <- function(organism=9606,max_pathway_size=80,...) {
  getBiocLiteLib('reactome.db')
  genesets <- list(...)
  path.list <- sapply(genesets,function(genelist) {
    list <- as.character(genelist)
    pathways <- subset(AnnotationDbi::select(reactome.db::reactome.db,keys=list,keytype="ENTREZID",columns=c("REACTOMEID","PATHNAME","ENTREZID")), ! is.na(REACTOMEID))
    pathways <- plyr::ddply(pathways,plyr::.(REACTOMEID), function(df) { if (reactome.total_count(df$REACTOMEID) <= max_pathway_size) { df } } )
    names(pathways) <- c('geneid','reactomeid','pathway')
    subset( merge(pathways,subset( convertEntrezIds(organism,genelist),uniprot %in% getUniprotIds(organism)),by='geneid') , select=c('reactomeid','pathway','geneid','genename','uniprot'))
  },USE.NAMES=T,simplify=F)
  for (name in names(path.list)) {
    attributes(path.list[[name]])$listname <- name
  }
  mergeddata <- Reduce(function(left,right) {
    rightlist <- subset( right, select=c('geneid','genename','uniprot','reactomeid'))
    if (is.null(attributes(left)$listname)) {
      names(rightlist) <- c(sapply( c('geneid','genename','uniprot'), function(x) { paste(x,attributes(right)$listname,sep='.')}),'reactomeid')
    }
    left <- merge(left,rightlist, by='reactomeid', suffixes=c(paste('.',attributes(left)$listname,sep=''),paste('.',attributes(right)$listname,sep='')))
    left
  }, path.list )
  mergeddata$db <- rep('reactomedb',nrow(mergeddata))
  rownames(mergeddata) <- mergeddata$reactomeid
  mergeddata$reactomeid<- NULL
  mergeddata
}