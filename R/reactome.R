reactome.total_count <- function(pathid) {
  nrow(AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID[pathid]))
}

#' @export
getCommonReactome <- function(organism=9606,max_pathway_size=80,...) {
  getBiocLiteLib('reactome.db')
  genesets <- list(...)
  path.list <- sapply(genesets,function(genelist) {
    list <- as.character(genelist)
    pathways <- subset(AnnotationDbi::select(reactome.db::reactome.db,keys=list,keytype="ENTREZID",columns=c("ENTREZID","REACTOMEID","PATHNAME")), ! is.na(REACTOMEID))
    pathways <- plyr::ddply(pathways,plyr::.(REACTOMEID), function(df) { if (reactome.total_count(df$REACTOMEID) <= max_pathway_size) { df } } )
    names(pathways) <- c('geneid','reactomeid','pathname')
    merge(pathways,subset( convertEntrezIds(organism,genelist),uniprot %in% getUniprotIds(organism)),by='geneid')
  },USE.NAMES=T,simplify=F)
  Reduce(function(left,right) {
    path.list[[left]] <- merge(path.list[[left]],subset( path.list[[right]], select=c('geneid','reactomeid','genename','uniprot')), by='reactomeid', suffixes=c(paste('.',left,sep=''),paste('.',right,sep='')))
    path.list[[left]]
  }, names(path.list) )
}