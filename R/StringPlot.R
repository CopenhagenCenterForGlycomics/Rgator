read_supplemental <- function() {
  uniprot.goa <- read.delim('gene_association.goa_human',header=F)[,c(2,5,7,9)]
  names(uniprot.goa)<- c('gene_id','go_id','Evidence','Ontology')
  #Actually uniprot ids
  uniprot.goa$Ontology <- mapvalues(uniprot.goa$Ontology,c('F','P','C'),c('MF','BP','CC'))
  mapped_ids_only <- getEntrezIds(9606,uniprot.goa$gene_id)
  mapped_ids <- data.frame(uniprot=names(mapped_ids_only),gene_id=as.character(mapped_ids_only))
  uniprot.goa$gene_id <- mapvalues(uniprot.goa$gene_id,mapped_ids$uniprot,mapped_ids$gene_id,warn_missing=F)
  return (uniprot.goa)
}

getLipidStringNetwork <- function() {
  supplemental=read_supplemental()
  rnaseq.t2enrichment <- getGOEnrichment(9606,c(),get_deseq_edge_ids(dge_analysis$t2,deseq_results$t2),ontology='BP',supplemental.terms=supplemental)
  diff.t2enrichment <- getGOEnrichment(9606,unique(subset(diff_data,t2_tcl <= -1 | t2_sec <= -1,select=c('uniprot','t2_tcl','t2_sec'))$uniprot),ontology='BP',supplemental.terms=supplemental)
  extra.rnaseq <- c('GO:0007009','GO:0015748','GO:0015866','GO:1902430','GO:0050435','GO:1901264')
  wanted_terms.rnaseq <- c(extra.rnaseq,'GO:0070328', 'GO:0071825', 'GO:0055088','GO:0006629')
  wanted_terms.diff <- c('GO:0042632','GO:0043691','GO:0070328','GO:0051006','GO:0010902','GO:0010984','GO:0071825','GO:0034375','GO:0034369','GO:0006638')
  wanted_ids <- unique( c('Q9Y5C1','P06727','Q86X29',getGOenrichmentGenes(rnaseq.t2enrichment,wanted_terms.rnaseq,9606)$uniprot,getGOenrichmentGenes(diff.t2enrichment,wanted_terms.diff,9606)$uniprot))
  wanted_ids <- intersect(wanted_ids,getUniprotIds(9606))
  network <- stringToNetwork(wanted_ids)
  return (list(ids=wanted_ids,network=network))
}