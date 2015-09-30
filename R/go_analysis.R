removeFactors <- function(df) {
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df
}

#' Get the genes that contributed to a particular enrichment of a set of terms
#'
#' @param   enrichment  Enrichment obtained from \code{\link{getGOEnrichment}}
#' @param   wanted_terms    vector of go ids that you wish to link to enriched genes
#' @param   organism    Organism that you wish to search for terms in
#' @return  Data frame with uniprot id, gene id and gene name
#' @export
getGOenrichmentGenes <- function(enrichment,wanted_terms=c(),organism=9606) {
  unique(convertEntrezIds(organism,unique(unlist(sapply(wanted_terms, function(go) {   geneIdsByCategory(enrichment)[[go]]  })))))
}

#' Get the GO terms associated with the given UniProt ids
#'
#' You can optionally supply a vector of GO terms to the method,
#' so that any children of the wanted terms will be considered as
#' the parent term. This is a cheap way of getting a slim GO set.
#' For example, if your wanted GO term was 'GO:0005829' (cytosol),
#' a gene annotated with the term 'GO:0044445' (cytosolic part)
#' will be mapped to 'GO:0005829' (cytosol).
#'
#' @param   organism  Organism to look for GO terms for
#' @param   uniprots  Vector of uniprot ids to retrieve GO terms for
#' @param   wanted    Top-level GO terms to search under
#' @param   ontology  GO ontology to search ('BP','MF','CC')
# @importFrom data.table rbindlist
# @importFrom AnnotationDbi Term
#' @export
getGOTerms <- function(organism=9606,uniprots,wanted=c(),ontology='BP') {
  all_terms <- retrieveGOTerms(organism,uniprots)
  if (length(wanted) > 0 && ontology %in% c('CC','BP','MF')) {
    go_ids <- sapply( wanted[wanted != 'GO:0005737'], function(x) { getGOChildren(x,ontology) }, USE.NAMES=T,simplify=F )
    if ( 'GO:0005737' %in% wanted ) {
      go_ids[['GO:0005737']] <- 'GO:0005737'
    }
    return ( do.call(rbind,sapply(names(go_ids),function(x) {
      expand.grid(go_id=x,term=AnnotationDbi::Term(x),
                  uniprot=unique(subset(all_terms, go_id %in% go_ids[[x]] & uniprot %in% uniprot )$uniprot))
    },simplify=F)) )
  } else {
    return(all_terms)
  }
}

# @importFrom GO.db GO.db
# @importFrom AnnotationDbi toTable
retrieveGOTerms <- function(organism=9606,uniprots) {
  getBiocLiteLib('GO.db')
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db')
  if (as.character(organism) == '559292') {
    organism <- 4932
  }
  if (as.character(organism) == '10029') {
    return (retrieveGOTerms.goa(organism,uniprots))
  }
  query_ids <- getEntrezIds(organism,uniprots)
  if (as.character(organism) == '4932') {
    entrez_ids <- query_ids
    entrez_mapping <- removeFactors(AnnotationDbi::toTable(revmap(org.Sc.sgd.db::org.Sc.sgdENTREZID)[query_ids]))
    query_ids <- unlist(entrez_mapping[1])
  }
  godb <- get( sub("\\.db","GO",organisms[as.character(organism)] ), envir=get(unlist( organisms[as.character(organism)])))
  terms <- removeFactors(AnnotationDbi::toTable(godb[query_ids]))
  names(terms) <- c('gene_id','go_id','Evidence','Ontology')
  if (as.character(organism) == '4932') {
    names(terms) <- c('systematic_name','go_id','Evidence','Ontology')
    terms <- merge(terms, data.frame(systematic_name=query_ids), by='systematic_name')
    terms <- merge(terms, entrez_mapping, by='systematic_name')
    terms <- uniqueframe( merge(terms, data.frame(uniprot=names(entrez_ids),gene_id=entrez_ids) , by='gene_id'))
    terms$uniprot <- tolower(terms$uniprot)
  } else {
    terms <- merge( terms, data.frame(gene_id=query_ids,uniprot=tolower(names(query_ids))), by='gene_id')
  }
  if (dim(terms)[1] > 0) {
    go_terms <- (select(GO.db::GO.db, as.character(terms$go_id),"TERM"))
    names(go_terms) <- c('go_id','term')
    terms <- merge( terms, go_terms, by='go_id')
    terms <- uniqueframe(terms)
    return (terms)
  }
  return (terms)
}

#' @rdname Rgator-deprecated
#' @export
GO_parents <- function(node = "GO:0008150", ontology = "BP") {
  .Deprecated('getGOParents',package='Rgator')
  getGOParents(node,ontology)
}


# @importFrom AnnotationDbi toTable
#' Get the parents of the given GO terms
#' @param node  GO node to examine
#' @param ontology  Ontology to look for parents in
#' @export
getGOParents <- function(node = "GO:0008150", ontology = "BP") {
    if (ontology == "BP") GOPARENTS <- GO.db::GOBPPARENTS
    if (ontology == "CC") GOPARENTS <- GO.db::GOCCPARENTS
    if (ontology == "MF") GOPARENTS <- GO.db::GOMFPARENTS
    kids <- node

    # initialize output
    out <- c(kids)
    # do the following until there are no more parents
    while (any(!is.na(kids))) {
        # Get the unique children of the parents (that aren't NA)
        #children <- unique(unlist(mget(parents[!is.na(parents)], envir=GOCHILDREN)))
        parents <- unique(unlist( removeFactors(AnnotationDbi::toTable(GOPARENTS[kids[!is.na(kids)]]))[2]))
        # append chldren to beginning of `out`
        # unique will keep the first instance of a duplicate
        # (i.e. the most recent child is kept)
        out <- unique(append(parents[!is.na(parents)], out))

        # children become the parents of the next generation
        kids <- parents
        if ('all' %in% kids) {
          kids <- c(NA)
        }
    }
    return(out)
}

#' @rdname Rgator-deprecated
#' @export
GO_children <- function(node = "GO:0008150", ontology = "BP") {
  .Deprecated('getGOChildren',package='Rgator')
  getGOChildren(node,ontology)
}

# @importFrom AnnotationDbi toTable
#' Get the parents of the given GO terms
#' @param node  GO node to examine
#' @param ontology  Ontology to look for children in
#' @export
getGOChildren <- function(node = "GO:0008150", ontology = "BP") {
    if (ontology == "BP") GOCHILDREN <- GO.db::GOBPCHILDREN
    if (ontology == "CC") GOCHILDREN <- GO.db::GOCCCHILDREN
    if (ontology == "MF") GOCHILDREN <- GO.db::GOMFCHILDREN
    parents <- node

    # initialize output
    out <- c(parents)
    # do the following until there are no more parents
    while (any(!is.na(parents))) {
        # Get the unique children of the parents (that aren't NA)
        #children <- unique(unlist(mget(parents[!is.na(parents)], envir=GOCHILDREN)))
        children <- unique(unlist(removeFactors(AnnotationDbi::toTable(GOCHILDREN[parents[!is.na(parents)]]))[1]))

        # append chldren to beginning of `out`
        # unique will keep the first instance of a duplicate
        # (i.e. the most recent child is kept)
        out <- unique(append(children[!is.na(children)], out))

        # children become the parents of the next generation
        parents <- children
    }
    return(out)
}

# @importFrom AnnotationDbi toTable
#' Perform an enrichment analysis on a given set of uniprot ids
#' @param organism  NCBI taxonomy identifier (e.g. 9606 for human)
#' @param uniprots  Vector of uniprot accessions to perform enrichment upon
#' @param query_ids Optional Entrez gene ids to use if there are no uniprot identifiers
#' @param universe  Vector of uniprot accessions to use as the universe. Defaults to all GO annotations for the organism
#' @param ontology  Ontology to perform enrichment in ('BP','CC','MF')
#' @param direction Direction to check for enrichment ('over' or 'under' for over- and under-representation respectively)
#' @param supplemental.terms  Flag to see if supplemental annotations should be retrieved
#'                            from UniprotGOA. Can also set to a dataframe containing the
#'                            terms.
#' @param conditional When set to TRUE, the enrichment will be conditional, taking into account relationships of terms
#' @return Enrichment object
#' @export
getGOEnrichment <- function(organism=9606,uniprots,query_ids=c(),universe=c(),ontology='BP',direction='over',supplemental.terms=F,conditional=TRUE) {
  getBiocLiteLib("GO.db")
  getBiocLiteLib("GOstats")
  getBiocLiteLib("GSEABase")
  organisms <- list('9606'='org.Hs.eg.db','10090'='org.Mm.eg.db','10116'='org.Rn.eg.db','7227'='org.Dm.eg.db','4932'='org.Sc.sgd.db','9823'='org.Ss.eg.db')
  dbname = NULL
  if ( as.character(organism) %in% names(organisms) ) {
    dbname = organisms[[as.character(organism)]]
    getBiocLiteLib(dbname)
    library(dbname,character.only=TRUE)
  }
  if (length(query_ids) < 1 & length(uniprots) > 0) {
    query_ids <- getEntrezIds(organism,uniprots)
  }
  if (is.null(dbname) || supplemental.terms == TRUE) {
    supplemental.terms <- downloadUniprotGOA(organism)
  }
  if (is.null(dbname) && !is.null(supplemental.terms)) {
    message("Don't have Gene IDs for this organism, using UniProt identifiers instead")
    query_ids <- uniprots
  }
  if (as.character(organism) == '4932' & length(query_ids) > 0) {
    entrez_ids <- query_ids
    entrez_mapping <- removeFactors(AnnotationDbi::toTable(revmap(org.Sc.sgdENTREZID)[query_ids]))
    query_ids <- unlist(entrez_mapping[1])
  }

  if (length(universe) < 1) {
    if (is.null(dbname)) {
      universe <- getUniprotIds(organism)
    } else {
      universe <- as.list( mappedkeys(get( sub("\\.db","GO",dbname ),asNamespace(dbname) )) )
    }
  } else {
    universe <- getEntrezIds(organism,universe)
  }
  library('GOstats')
  if (is.na(supplemental.terms) || supplemental.terms == FALSE) {
  params <- new('GOHyperGParams',
              geneIds=query_ids,
              universeGeneIds=unlist(universe),
              ontology=ontology,
              pvalueCutoff=0.05,
              conditional=conditional,
              testDirection=direction,
              annotation=as.character(organisms[as.character(organism)])
             )
  } else {
    library("GSEABase")
    target_db = supplemental.terms
    column_name = "uniprot"
    if (! is.null(dbname) ) {
      target_db = rbind( target_db, removeFactors( as.data.frame(AnnotationDbi::toTable( get( sub("\\.db","GO", dbname),asNamespace(dbname) ) ))) )
      column_name = "gene_id"
    }
    super_terms = unique(target_db[target_db$Ontology == ontology,])
    frame=GOFrame(data.frame(super_terms$go_id, super_terms$Evidence, super_terms[[column_name]]),organism=as.character(organism))
    allFrame=GOAllFrame(frame)
    gsc <- GeneSetCollection(allFrame, setType = GOCollection())
    params <- GSEAGOHyperGParams(name="Supplemental GO annotation",
                             geneSetCollection=gsc, geneIds = query_ids, universeGeneIds = unlist(universe),
                             ontology = ontology, pvalueCutoff = 0.05, conditional = conditional, testDirection =  direction)
  }

  # Reference for info on multiple testing correction, and why we don't do it:
  # https://stat.ethz.ch/pipermail/bioconductor/2008-January/020690.html

  hgOver <- hyperGTest(params)
  return (hgOver)
}

retrieveGOTerms.goa <- function(organism,uniprots) {
  goa <- downloadUniprotGOA(organism)
  goa[goa$uniprot %in% toupper(uniprots),]
}

#' Download annotations from the UniprotGOA site
#' @param   organism  NCBI taxonomy id (e.g. 9606) to get terms for
#' @return  Data frame with gene id, go_id, Evidence and Ontology columns
#' @export
downloadUniprotGOA <- function(organism=9606) {
  organism <- as.character(organism)
  organisms <- list('9606'='human','10090'='mouse','10116'='rat','7227'='fly','4932'='yeast')
  proteomes <- list('10029'='264824.C_griseus', '9823'='35497.S_scrofa.goa')
  if (organism %in% names(organisms)) {
    species<-organisms[[organism]]
    uniprot.goa <- cacheFile( paste("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/",toupper(species),"/gene_association.goa_",species,".gz",sep=''),paste("goa-",species,sep=""),gzip=T,stringsAsFactors=F,comment.char='!')[,c(2,5,7,9)]
  }
  if (organism %in% names(proteomes)) {
    species<-proteomes[[organism]]
    uniprot.goa <- cacheFile( paste("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/",species,".goa",sep=''),paste("goa-",species,sep=""),gzip=F,stringsAsFactors=F,comment.char='!')[,c(2,5,7,9)]
  }
  if (is.null(uniprot.goa)) {
    message("Invalid organism")
    return ()
  }
  names(uniprot.goa)<- c('gene_id','go_id','Evidence','Ontology')
  #This frame is actually indexed by uniprot ids, but call the column gene_id anyway to avoid renaming
  uniprot.goa$Ontology <- plyr::mapvalues(uniprot.goa$Ontology,c('F','P','C'),c('MF','BP','CC'))
  mapped_ids_only <- getEntrezIds(organism,uniprot.goa$gene_id)
  if (is.null(mapped_ids_only)) {
    names(uniprot.goa) <- c('uniprot','go_id','Evidence','Ontology')
    return (uniprot.goa)
  }
  mapped_ids <- data.frame(uniprot=names(mapped_ids_only),gene_id=as.character(mapped_ids_only))
  uniprot.goa$gene_id <- plyr::mapvalues(uniprot.goa$gene_id,mapped_ids$uniprot,mapped_ids$gene_id,warn_missing=F)
  return (uniprot.goa)
}

#' Get a list of potentially cytosolic proteins from uniprot identifiers
#' Simple check that sees if there are only cytosolic related terms associated
#' with a given identifier
#' @param   organism  NCBI taxonomy id (e.g. 9606) to get terms for
#' @param   uniprots  UniProt identifiers to check for cytosolic
#' @return  Vector of UniProt identifiers
#' @export
getCytosolic <- function(organism,uniprots) {
  tabled_terms <- getGOTerms(organism,unique(uniprots),wanted=c('GO:0030054','GO:0005886','GO:0012505','GO:0005783','GO:0005794','GO:0005764','GO:0016020','GO:0005739','GO:0005634','GO:0005576','GO:0005773','GO:0005829','GO:0005856','GO:0005737'),ontology='CC')
  potential_cytosol <- as.character(unique(subset(tabled_terms,term %in% c('nucleus','vacuole','cytosol','cytoskeleton','cytoplasm'))$uniprot))
  potential_extracellular <- subset(tabled_terms, ! term %in% c('nucleus','vacuole','cytosol','cytoskeleton','cytoplasm'))$uniprot
  potential_cytosol[ ! potential_cytosol %in% potential_extracellular ]
}

#' @rdname Rgator-deprecated
#' @export
tableGOTerms <- function(organism,uniprot,wanted=c('GO:0030054','GO:0005886','GO:0012505','GO:0005783','GO:0005794','GO:0005764','GO:0016020','GO:0005739','GO:0005634','GO:0005576','GO:0005773','GO:0005829','GO:0005856','GO:0005737'),ontology='CC') {
  .Deprecated('getGOTerms',package='Rgator')
  getGOTerms(organism,uniprot,wanted,ontology)
}