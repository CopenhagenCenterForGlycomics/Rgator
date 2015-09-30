#' Replace CDD identifiers with cluster identifiers
#'
#' By default, methods on the CDD data sources use the
#' cddid as the identifier for a domain. We can use the
#' cluster identifier as an alternative if we wish to
#' group similar domains.
#'
#' This method connects to the CDD to download the
#' family to superfamily links for the database, and
#' caches this file locally.
#'
#'  @param   cddids  Vector of cddids to look up
#'  @return  vector of identifiers
#'  @export
cddidToSuperfamily <- function(cddids) {
  cdd_superfamily_links <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links","cdd-family-superfamily")[,c(1,3)]
  names(cdd_superfamily_links) <- c('cddid','clusterid')
  for(i in 1:nrow(cdd_superfamily_links)) {
    cddids[ cddids==cdd_superfamily_links[[i,'cddid' ]] ] <- cdd_superfamily_links[[ i , 'clusterid'  ]]
  }
  cddids
}

#' Retrieve the names for a vector of cdddids
#'
#' Names are retrieved from the CDD database for all domains
#' and a data frame with columns \code{dom}, \code{id}, \code{short}, \code{long} is returned
#'  @param   cddids   Vector of cddids to lookup names for
#'  @param   data.frame Domain names
#'  @export
getCddNames <- function(cddids) {
  cddid_all <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz","cddid-all",gzip=T)
  names(cddid_all) <- c('id','dom','short','long')
  merge(data.frame(dom=cddids),cddid_all,by='dom',all.x=T)
}

#' Get domain sets - Deprecated
#' @rdname Rgator-deprecated
#' @export
getDomainSets <- function( inputsites, sitecol, domaindata,  max_dom_proportion=0.81, stem_distance=100  ) {
  .Deprecated('calculateDomainSets',package='Rgator')
  calculateDomainSets(inputsites, sitecol, domaindata,  max_dom_proportion, stem_distance)
}

#' @export
downloadDomains <- function(...) {
  organism = as.character(list(...))
  if ("9606" %in% organism) {
    set = downloadDataset('http://glycodomain-data.glycocode.com/data/latest/fulldomains/',list(type='gatorURL',title='fulldomains'))
    set
  }
  sapply( organism[organism != '9606'], function(org) {
    set = downloadDataset(paste('http://glycodomain-data.glycocode.com/data/latest/domains.',org,'/',sep=''),list(type='gatorURL',title=paste('domains.',org,sep='')))
    set
  })
}

#' @export
getGOterms.interpro <- function(interpro=NULL,wanted=c(),ontology='BP') {
  interpro_go <- cacheFile('ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go','interpro-go',comment.char='!')$V1
  interpro_ids <- trimws(substr(interpro_go,10,18))
  go_ids <- trimws(substr(interpro_go,nchar(interpro_go)-10,nchar(interpro_go)))
  go_mapping <- data.frame(interpro=interpro_ids,go_id=go_ids)
  ontologies <- AnnotationDbi::select(GO.db::GO.db, keys=trimws(unique(go_ids)),keytype='GOID',c('ONTOLOGY','GOID'))
  names(ontologies) <- c('go_id','Ontology')
  go_mapping <- merge(go_mapping,ontologies,by='go_id',all.x=T)
  if ( !is.null(interpro) ) {

    if (length(wanted) > 0 && ontology %in% c('CC','BP','MF')) {
      go_ids <- sapply( wanted, function(x) { getGOChildren(x,ontology) }, USE.NAMES=T,simplify=F )
      return ( as.data.frame(data.table::rbindlist(sapply(names(go_ids),function(x) {
        expand.grid(go_id=x,term=AnnotationDbi::Term(x),
                    interpro=unique(go_mapping[ go_mapping$go_id %in% go_ids[[x]] & go_mapping$interpro %in% interpro,'interpro' ] ))
      },simplify=F)) ))
    } else {
      return (go_mapping[go_mapping$interpro %in% interpro,])
    }


  } else {
    go_mapping
  }
}

#' @export
downloadInterproDomains <- function(...) {
  organism = as.character(list(...))
  interpro_classes <- cacheFile('ftp://ftp.ebi.ac.uk/pub/databases/interpro/53.0/entry.list','interpro-classes')
  starts = c( sort(sapply(c('Active_site','Binding_site','Conserved_site','Domain','Family','PTM','Repeat'), function(clazz) { which(interpro_classes$V1 == clazz) })), nrow(interpro_classes) )
  classes <- list()
  for (row in 1:(length(starts)-1)) {
    class_data = data.frame(V1=interpro_classes[ (starts[row]+1):(starts[row+1]-1), ])
    class_data$dom = substr(class_data$V1,1,9)
    class_data$description <- substr(class_data$V1,11,nchar(class_data$V1))
    class_data$class <- names(starts)[row]
    classes[[ names(starts)[row] ]] <- class_data
  }
  annotations <- do.call(rbind,classes)
  sapply( organism, function(org) {
    setname = paste('interpro.',org,sep='')
    downloadDataset(paste('http://glycodomain-data.glycocode.com/data/latest/interpro.53.0.',org,'/',sep=''),list(type='gatorURL',title=setname))
    data.env = getDataEnvironment()
    data.env[[ paste(setname,'annotated',sep='.') ]] <- merge(get(setname),annotations,all.x=T,by='dom')[,c('uniprot','dom','start','end','class','description')]
    saveEnvironment()
  })
}

downloadDisulfides <- function(...) {
  organism = as.character(list(...))
  for (org in organism) {
    disulfides <- cacheUniprotFile(paste("http://www.uniprot.org/uniprot/?query=annotation%3A(type%3Adisulfid)&format=tab&columns=id,feature(DISULFIDE%20BOND)&fil=organism%3A",org,sep=''),paste('disulfide-uniprot-',org,sep=''),header=T);
    extracted = plyr::llply( strsplit(disulfides$Disulfide.bond,'[0-9\\.]; '), function(el) {
      res = Map( function(splt) {
        if ("Interchain" %in% splt) {
          c( as.numeric(splt[c(2,3)]), "interchain")
        } else {
          c( as.numeric(splt[c(2,3)]), "intrachain")
        }
      }, strsplit(el,' '))
      return (res)
    })
    names(extracted) <- tolower(disulfides$Entry)
    fixed = plyr::ldply(extracted,function(prot) {
      res = as.data.frame(matrix(unlist(prot),ncol = 3,byrow=T))
    })
    names(fixed) <- c('uniprot','start','end','type')
    attributes(fixed)$version <- paste('UniProt',attributes(disulfides)$version,sep='_')
    data.env = getDataEnvironment()
    data.env[[ paste('disulfide.uniprot.',org,sep='') ]] <- fixed
    saveEnvironment()
    updateDataVersions()
  }
}

downloadProcessing <- function(...) {
  organism = as.character(list(...))
  for (org in organism) {
    processing <- cacheUniprotFile(paste("http://www.uniprot.org/uniprot/?query=annotation%3A(type%3Apeptide)+OR+annotation%3A(type%3Apropep)&format=tab&columns=id,feature(PROPEPTIDE),feature(PEPTIDE)&fil=organism%3A",org,sep=''),paste('processing-uniprot-',org,sep=''),header=T);
    extracted.propep = plyr::llply( strsplit(processing$Propeptide,'PROPEP '), function(el) {
      res = Map( function(splt) {
        c( as.numeric(splt[c(1,2)]), "propep" )
      }, strsplit(el,' '))
      return (res)
    })
    names(extracted.propep) <- tolower(processing$Entry)
    extracted.pep = plyr::llply( strsplit(processing$Peptide,'PEPTIDE '), function(el) {
      res = Map( function(splt) {
        c( as.numeric(splt[c(1,2)]), "peptide" )
      }, strsplit(el,' '))
      return (res)
    })
    names(extracted.pep) <- tolower(processing$Entry)
    fixed.pep = plyr::ldply(Filter(length,extracted.pep),function(prot) {
      res = as.data.frame(matrix(unlist(prot),ncol = 3,byrow=T))
    })
    names(fixed.pep) <- c('uniprot','start','end','type')
    fixed.propep = plyr::ldply(Filter(length,extracted.propep),function(prot) {
      res = as.data.frame(matrix(unlist(prot),ncol = 3,byrow=T))
    })
    names(fixed.propep) <- c('uniprot','start','end','type')
    fixed = rbind(fixed.pep,fixed.propep)
    fixed = fixed[!is.na(fixed$start),]
    attributes(fixed)$version <- paste('UniProt',attributes(processing)$version,sep='_')
    data.env = getDataEnvironment()
    data.env[[ paste('processing.uniprot.',org,sep='') ]] <- fixed
    saveEnvironment()
    updateDataVersions()
  }
  ;
}

#' Divide up sets of sites based on their site context
#'
#' It is useful to be able to classify sites so that you can see where the sites
#' occur in relation to other features on a protein. One approach is to look at
#' the conserved domains (often a proxy for the degree of amino acid conservation
#' in a given region), and then classify sites based upon the location they are found
#' in relation to these conserved domains.
#'
#' Sequences are retrieved from UniProt, and a domain data frame \code{\link{downloadDomains}} is required
#' Some conserved domains span the entire length of the protein, and as such they may not
#' be informative. To filter these domains, the \code{max_dom_proportion} parameter is used
#' to filter the domain data so that they are not used in the classification process (\code{real} domains).
#'
#' \strong{Basic region classification}
#' Sites that occur within the domain (\code{inside} domains) are chosen based on the site position
#' being within the start and end positions for the domain, and are not within the predicted transmembrane
#' region. All combinations of sites and domains are given.
#'
#' Sites that occur outside the domain (\code{outside} domains) are all the remaining sites that are not
#' classified as \code{inside}. All combinations of sites and domains in the protein are given. The domains
#' that are just outside in this list aren't particularly useful as they are largely a recaptiulation of the
#' domains that exist within the protein.
#'
#' Sites that occur both N- and C- terminal from a domain and are \code{outside}, are defined as being \code{between}
#' domains. All sets of N- and C- terminal domains for a given site are given. For example, if a site is
#' located between two domains (imagine a linker domain), then there are domains both N- and C- terminal of the site.
#' If the site, however, is located at the C-terminal, then it will only have a list of N- domains.
#'
#' \strong{Stem regions}
#' Type I stem regions (\code{stem.typei}) are single-pass TM proteins having a transmembrane region
#' at the C-terminal of the protein. Sometimes conserved domains are defined on
#' top of transmembrane regions, so to deal with this, Type I transmembrane
#' regions are defined as regions where the closest (or next-closest) C-terminal
#' domain to the site is a transmembrane region. Only sites lying within
#' \code{stem_distance} amino acids are returned.
#'
#'
#' Type II transmembrane proteins are defined as being single-pass TM proteins having a transmembrane region at the N-terminal
#' of the protein. Given that transmembrane regions may be broken down into smaller transmembrane
#' regions by the prediction algorithms, the definition of a Type II is given as a set of N-terminal
#' transmembranes that overlap a signlap peptide prediction, or a single transmembrane region on a
#' protein (without any signal peptide).
#'
#' Type II stem regions (\code{stem.typeii})are defined as \code{between} sites that occur within \code{stem_distance} amino acids of the
#' transmembrane region of the protein. All domains that are within \code{stem_distance} amino acids of the site,
#' as well as the closest transmembrane region or signal peptide are returned.
#'
#' SIGNALP stem regions (\code{stem.signalp}) aren't strictly stem regions, but are defined as the
#' region between the Signal peptide cleavage site and any domain. Although the
#' signal peptide is usually cleaved, there are cases where the signal peptide
#' does not get cleaved. They are defined in the same way as Type II stem regions,
#' except there N-terminal transmembrane region is instead a signal peptide.
#'
#' \strong{Other region types}
#' Interdomain regions (\code{interdomain}, also described as linker regions)
#' contain sites that are part of the \code{between} set (i.e. between two domains)
#' but do not fit into the categories of stem regions.
#'
#' N- or C-terminal regions (\code{norc}), contain sites that are outside, but
#' do not fit into any of the other categories of site.
#'
#' @param   inputsites          Data frame with a uniprot column and a column containing sites
#' @param   sitecol             Column to look for site data in
#' @param   domaindata          Data frame containing all the domain data
#' @param   max_dom_proportion  Maximum proportion of total protein length domain can take
#' @param   stem_distance       Maximum distance from the transmembrane region before a region is no longer a stem.
#' @return  named list of data.frames containing sites from each of the regions (\code{inside}, \code{outside}, \code{stem},\code{stem.typei},\code{stem.typeii},\code{stem.signalp},\code{interdomain},\code{norc})
#' @seealso \code{\link{downloadDomains}}
#'
#'
# @importFrom plyr .
# @importFrom plyr ddply
#' @export
calculateDomainSets <- function( inputsites, sitecol, domaindata, max_dom_proportion=0.81, stem_distance=100  ) {
  message("Retrieving Uniprot sequences")
  seqdat <- getUniprotSequences(unique(inputsites$uniprot))
  message("Retrieved sequences")
  domdat <- merge(merge(domaindata,subset(inputsites,! is.na(inputsites[[sitecol]])),by='uniprot',allow.cartesian=TRUE),seqdat,by='uniprot',allow.cartesian=TRUE)
  domdat$aalength <- nchar(domdat$sequence)
  domdat$sitekey <- paste(domdat$uniprot,'-',domdat[[sitecol]],sep='')
  domdat$sequence <- NULL
  real <- subset( domdat, ((as.numeric(end) - as.numeric(start)) / aalength < max_dom_proportion ))
  inside <- subset( subset (  real ,  ( (as.numeric(real[[sitecol]]) >= as.numeric(start)) & (as.numeric(real[[sitecol]]) <= as.numeric(end))  )  ), ! grepl("tmhmm",dom))
  outside <- subset ( real , ! sitekey %in% inside$sitekey & dom != 'tmhmm-outside' & dom != 'tmhmm-inside' )
  sitekeys_nterm <- unique(subset( outside, as.numeric(outside[[sitecol]]) < as.numeric(start) )$sitekey)
  sitekeys_cterm <- unique(subset( outside, as.numeric(outside[[sitecol]]) > as.numeric(end) )$sitekey)
  between <- subset( outside, sitekey %in% intersect(sitekeys_nterm, sitekeys_cterm ) )
  message("Identifying type II proteins")
  typeii <- unique(plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- unique(subset(input,select=c('dom','start','end','sitekey')))
    signals <- subset(df, dom == 'SIGNALP')
    tms <- subset(df, dom == "tmhmm-TMhelix")
    if (dim(signals)[1] < 1) {
      if (dim(tms)[1] == 1) {
        return (input)
      }
      return ()
    }
    signal_start <- signals$start[1]
    signal_end <- signals$end[1]
    tms <- subset(df, dom == "tmhmm-TMhelix" & as.numeric(start) < as.numeric(signal_end) )
    other_tms <- subset(df, dom == "tmhmm-TMhelix" & as.numeric(start) > as.numeric(signal_end) )

    if (dim(tms)[1] > 0 & dim(other_tms)[1] == 0) {
      return (input)
    }
    return ()
  },.progress="text")$uniprot)
  message("Identifying type II stems")
  stem_typeii <- plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the N-terminal domains, our site is C-terminal of the domain
    filtered <- subset(df,siteend>0)
    filtered <- filtered[order(filtered$siteend),]
    # We have a SIGNALP or TMHELIX -- something...
    # This is a type II transmembrane
    if (( (filtered$dom[1] == "SIGNALP") | grepl("tmhmm-TMhelix",filtered$dom[1]) ) & (filtered$uniprot[1] %in% typeii )) {
      # Grab all the C-terminal domains within the stem distance of the site
      # Grab all the N-terminal domains within the stem distance of the site, and closer than the closest SIGNALP or TM (i.e the closest SIGNALP or TM)
      wanted <- subset(df, ((startsite > 0 & startsite <= stem_distance) | (siteend > 0 & siteend <= stem_distance & siteend <= filtered$siteend[1] )))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
  },.progress="text")
  message("Identifying Signalp stems")
  signalp_stem <- plyr::ddply(between,plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the N-terminal domains, our site is C-terminal of the domain
    filtered <- subset(subset(df,siteend>0),dom != 'tmhmm-TMhelix')
    filtered <- filtered[order(filtered$siteend),]
    # We have a SIGNALP -- something...
    # If we've got a signalp then it is secreted


    if ( dim(filtered)[1] > 0 & ( (filtered$dom[1] == "SIGNALP") ) & (! filtered$uniprot[1] %in% typeii )) {
      # Anything with a transmembrane helix is not secreted
      wanted <- subset(df, ((siteend > 0 & siteend <= stem_distance & siteend <= filtered$siteend[1] )))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
  },.progress="text")
  message("Identifying type I stems")
  stem_typei <- plyr::ddply(subset(between, ! sitekey %in% signalp_stem$sitekey),plyr::.(sitekey),function(input) {
    df <- input
    df$start <- as.numeric(df$start)
    df$siteend <- as.numeric(df[[sitecol]]) - as.numeric(df$end)
    df$startsite <- as.numeric(df$start) - as.numeric(df[[sitecol]])

    # All the C-terminal domains, our site is N-terminal of the domain
    # This is a type I transmembrane
    filtered <- subset(df,startsite>0)
    filtered <- filtered[order(filtered$startsite),]
    # We have a something -- TMHELIX

    # If we have a domain spanning the transmembrane, we want to call this a type I transmembrane too

    if ( ((dim(filtered)[1] >= 1) & ( filtered$dom[1] == "tmhmm-TMhelix" )) | ((dim(filtered)[1] >= 2) & (filtered$dom[2] == "tmhmm-TMhelix") & (as.numeric(filtered$end[1]) > as.numeric(filtered$end[2]) ) )) {
      wanted <- subset(df, ((startsite > 0 & startsite <= stem_distance & startsite <= filtered$startsite[1] ) | (siteend > 0 & siteend <= stem_distance)))
      wanted$siteend <- NULL
      wanted$startsite <- NULL
      return (wanted)
    }
    return ()
  },.progress="text")

  interdomain <- subset(between, ! sitekey %in% stem_typei$sitekey & ! sitekey %in% stem_typeii$sitekey & ! sitekey %in% signalp_stem$sitekey )
  norc <- subset(outside, ! sitekey %in% between$sitekey )
  #Stem = Betweeen where closest N-terminal = SIGNALP/TMHMM
  #ddply between by sitekey if (site - end), sort asc [1] $dom == tmhmmm/signalp return df
  #                         if (start - site), sort asc [1] $dom == tmhmm/signalp return df
  #                         else return empty
  return ( list( all=domdat, real=real, inside=inside, outside=outside, between=between, stem=rbind(stem_typei,stem_typeii,signalp_stem), stem.typeii=stem_typeii, stem.typei=stem_typei, stem.signalp=signalp_stem, interdomain=interdomain, norc=norc  )  )
}