#' Plot a String network from the given seed UniProt accession identifiers
#'
#' At a given evidence threshold, we can retrieve the String database relations, and add them
#' to a graph object, which is then rendered
#'  @param  organism  NCBI taxonomy ID for which to retrieve the String network
#'  @param  accessions  Vector of UniProt accessions to retrieve the String network for
#'  @param  get.neighbours  Boolean flag (TRUE / FALSE) indicating if immediate neighbours of the
#'                          given accessions should also be added to the graph
#'  @param  threshold   Minimum threshold for the combined score criteria for addition to the graph
#'  @return List containing two elements, plot: Containing the plot, and coords: Containing the co-ordinates of each protein
#'  @export
generateStringNetwork <- function(organism=9606,accessions=c('q14118','P24043'),get.neighbours=F,threshold=700) {
  ggnet(getStringNetwork(organism,accessions,get.neighbours,threshold))
}

# @importFrom STRINGdb STRINGdb
# @importFrom AnnotationDbi select
# @importFrom intergraph asDF
# @importFrom plyr mapvalues
# @importFrom network network set.edge.attribute set.vertex.attribute network.vertex.names
#' Get a String network from the given seed UniProt accession identifiers
#'
#' At a given evidence threshold, we can retrieve the String database relations, and add them
#' to a graph object. This graph can then be drawn using the \code{\link{ggnet}} method
#'  @param  organism  NCBI taxonomy ID for which to retrieve the String network
#'  @param  accessions  Vector of UniProt accessions to retrieve the String network for
#'  @param  get.neighbours  Boolean flag (TRUE / FALSE) indicating if immediate neighbours of the
#'                          given accessions should also be added to the graph
#'  @param  threshold   Minimum threshold for the combined score criteria for addition to the graph
#'  @return A network class object
#'  @export
#'  @seealso \code{\link[network]{network}} \code{\link[STRINGdb]{STRINGdb}} \code{\link{ggnet}}
getStringNetwork <- function(organism=9606,accessions=c(),get.neighbours=F,threshold=700) {
  getBiocLiteLib('STRINGdb')
  require('igraph')
  string_db <- STRINGdb::STRINGdb$new(species=organism,version="9_05",score_threshold=threshold)

  id_mappings <- prepare_string_ids(string_db,organism,accessions,get.neighbours)

  subnetwork <- intergraph::asDF(string_db$get_subnetwork(id_mappings$STRING_id))

  subnetwork$vertexes$name <- plyr::mapvalues(subnetwork$vertexes$name,id_mappings$STRING_id,tolower(id_mappings$uniprot),warn_missing=F)
  subnetwork$edges$V1 <- plyr::mapvalues(subnetwork$edges$V1, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  subnetwork$edges$V2 <- plyr::mapvalues(subnetwork$edges$V2, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  net <- network::network(subset(subnetwork$edges,select=c(1,2),!is.na(V1) & !is.na(V2) ),directed=F)

  network::set.edge.attribute(net,attrname="combined",value=subnetwork$edges$combined_score)
  network::set.edge.attribute(net,attrname="experimental",value=subnetwork$edges$experimental)
  network::set.edge.attribute(net,attrname="type",value=rep(NA,nrow(subnetwork$edges)))

  # Sure, a "rep" call here will do the job, but we might want to change the class in the future
  network::set.vertex.attribute(net,attrname='class',value=as.character(sapply(network::network.vertex.names(net),function(id) { "Protein" })))
  network::set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network::network.vertex.names(net),function(id) { toupper(id) })))
  return (net)
}

prepare_string_ids <- function(string_db,organism=9606,accessions=c(),get.neighbours=F) {
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  dbname<-organisms[[as.character(organism)]]
  getBiocLiteLib(dbname)
  id_mappings <- AnnotationDbi::select(get(dbname),columns=c('UNIPROT','ENSEMBLPROT'),keys=toupper(accessions),keytype='UNIPROT')
  id_mappings$ENSEMBLPROT <- sapply ( id_mappings$ENSEMBLPROT, function(x) { paste(organism,'.',x,sep='') })
  names(id_mappings) <- c('uniprot','STRING_id')
  all_string_proteins <- string_db$get_proteins()$protein_external_id

  # The id_mappings frame determine which IDs we will be retrieving from the String database.

  id_mappings <- unique(subset(id_mappings, STRING_id %in% all_string_proteins ))

  # If we wish to get the neighbours for a set of IDs, we use the get_neighbours method from the STRINGdb
  # object, retrieve out the ENSEMBL identifiers, and then convert them across to UniProt ids, appending
  # to the id_mappings data frame.
  if (get.neighbours) {
    in_subnetwork <- intersect( id_mappings$STRING_id, intergraph::asDF(string_db$get_subnetwork(id_mappings$STRING_id))$vertexes$name)


    neighbours <- string_db$get_neighbors(in_subnetwork)
    neighbour_mappings <- AnnotationDbi::select(get(dbname),columns=c('UNIPROT','ENSEMBLPROT'),keys=gsub(paste(organism,'.',sep=''),"",neighbours),keytype='ENSEMBLPROT')
    neighbour_mappings$ENSEMBLPROT <- sapply ( neighbour_mappings$ENSEMBLPROT, function(x) { paste(paste(organism,'.',sep=''),x,sep='') })
    names(neighbour_mappings) <- c('STRING_id','uniprot')
    neighbour_mappings <- subset(neighbour_mappings,select=c('uniprot','STRING_id'))
    id_mappings <- rbind(id_mappings,neighbour_mappings)
  }
  return(id_mappings)
}

getStringClusters <- function(organism=9606,accessions=c(),version="9_1",get.neighbours=F,threshold=700,total.clusters=NA,algorithm=NA,...) {
  getBiocLiteLib('STRINGdb')
  require('igraph')
  string_db <- STRINGdb::STRINGdb$new(species=organism,version=version,score_threshold=threshold)

  id_mappings <- prepare_string_ids(string_db,organism,accessions,get.neighbours)

  if (is.na(algorithm)) {
    algorithm <- igraph::fastgreedy.community
  }
  community.data <- algorithm( string_db$get_subnetwork(id_mappings$STRING_id), ... )
  if (! is.na(total.clusters) && igraph::is.hierarchical(community.data) ) {
    memb = cutat(community.data,total.clusters)
  } else {
    memb = membership(community.data)
  }
  clusters = NULL
  for (i in 1:max(memb)) {
      clusters[[i]] = names(membership(community.data)[membership(community.data) == i])
  }
  clusters <- lapply(clusters, function(clust) { plyr::mapvalues(clust,id_mappings$STRING_id,tolower(id_mappings$uniprot),warn_missing=F) } )
  return (clusters)
}