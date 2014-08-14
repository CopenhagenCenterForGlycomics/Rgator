# @importFrom STRINGdb STRINGdb
# @importFrom AnnotationDbi select
# @importFrom intergraph asDF
# @importFrom plyr mapvalues
# @importFrom network network set.edge.attribute set.vertex.attribute network.vertex.names
getStringNetwork <- function(organism=9606,uniprots=c('q14118','P24043'),get.neighbours=F,threshold=700) {
  getBiocLiteLib('STRINGdb')
  string_db <- STRINGdb::STRINGdb$new(species=organism,version="9_05",score_threshold=threshold)
  all_string_proteins <- string_db$get_proteins()$protein_external_id
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  id_mappings <- AnnotationDbi::select(get(dbname),columns=c('UNIPROT','ENSEMBLPROT'),keys=toupper(uniprots),keytype='UNIPROT')
  id_mappings$ENSEMBLPROT <- sapply ( id_mappings$ENSEMBLPROT, function(x) { paste(organism,'.',x,sep='') })
  names(id_mappings) <- c('uniprot','STRING_id')
  id_mappings <- unique(subset(id_mappings, STRING_id %in% all_string_proteins ))
  if (get.neighbours) {
    neighbours <- string_db$get_neighbors(id_mappings$STRING_id)
    neighbour_mappings <- AnnotationDbi::select(get(dbname),columns=c('UNIPROT','ENSEMBLPROT'),keys=gsub(paste(organism,'.',sep=''),"",neighbours),keytype='ENSEMBLPROT')
    neighbour_mappings$ENSEMBLPROT <- sapply ( neighbour_mappings$ENSEMBLPROT, function(x) { paste(paste(organism,'.',sep=''),x,sep='') })
    names(neighbour_mappings) <- c('STRING_id','uniprot')
    neighbour_mappings <- subset(neighbour_mappings,select=c('uniprot','STRING_id'))
    id_mappings <- rbind(id_mappings,neighbour_mappings)
  }
  subnetwork <- intergraph::asDF(string_db$get_subnetwork(id_mappings$STRING_id))

  subnetwork$vertexes$name <- plyr::mapvalues(subnetwork$vertexes$name,id_mappings$STRING_id,tolower(id_mappings$uniprot),warn_missing=F)
  subnetwork$edges$V1 <- plyr::mapvalues(subnetwork$edges$V1, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  subnetwork$edges$V2 <- plyr::mapvalues(subnetwork$edges$V2, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  net <- network::network(subset(subnetwork$edges,select=c(1,2),!is.na(V1) & !is.na(V2) ),directed=F)
  network::set.edge.attribute(net,attrname="combined",value=subnetwork$edges$combined_score)
  network::set.edge.attribute(net,attrname="experimental",value=subnetwork$edges$experimental)
  network::set.edge.attribute(net,attrname="type",value=rep(NA,nrow(subnetwork$edges)))
  network::set.vertex.attribute(net,attrname='class',value=as.character(sapply(network::network.vertex.names(net),function(id) { "Protein" })))
  network::set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network::network.vertex.names(net),function(id) { toupper(id) })))
  return (net)
}