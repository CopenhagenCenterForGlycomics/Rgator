#' Given sets of genes, find the common pathways that these genes belong to
#'
#' This looks within ReactomeDB, and Graphite (encompassing many other databases)
#' @param organism          NCBI taxonomy ID for the organism that is being studied
#' @param max_pathway_size  Maximum size of the pathway that should be checked
#' @param ...               Lists of Entrez gene identifiers
#' @return Data frame with columns pathway, db (indicating source DB) and the gene identifiers that are common to that pathway
#' @examples
#' findCommonPathways(organism=9606,max_pathway_size=1000,seta=c('1605'),setb=c('1605'))
#' @export
findCommonPathways <- function(organism=9606,max_pathway_size=80,...) {
  rbind( findCommonPathways.reactome(organism,max_pathway_size,...), findCommonPathways.graphite(organism,max_pathway_size,...), use.names=T )
}

reactome.total_count <- function(pathid) {
  nrow(AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID[pathid]))
}

#' Given sets of genes, find the common pathways that these genes belong to within ReactomeDB
#'
#' @param organism          NCBI taxonomy ID for the organism that is being studied
#' @param max_pathway_size  Maximum size of the pathway that should be checked
#' @param ...               Lists of Entrez gene identifiers
#' @return Data frame with columns pathway, db (indicating source DB) and the gene identifiers that are common to that pathway
#' @examples
#' findCommonPathways(organism=9606,max_pathway_size=1000,seta=c('1605'),setb=c('1605'))
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