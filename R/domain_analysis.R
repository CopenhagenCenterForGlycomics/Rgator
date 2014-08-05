#' @export
cddidToSuperfamily <- function(cddids) {
  cdd_superfamily_links <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links","cdd-family-superfamily")[,c(1,3)]
  names(cdd_superfamily_links) <- c('cddid','clusterid')
  for(i in 1:nrow(cdd_superfamily_links)) {
    cddids[ cddids==cdd_superfamily_links[[i,'cddid' ]] ] <- cdd_superfamily_links[[ i , 'clusterid'  ]]
  }
  cddids
}

#' @export
getCddNames <- function(cddids) {
  cddid_all <- cacheFile("ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz","cddid-all",gzip=T)
  names(cddid_all) <- c('id','dom','short','long')
  merge(data.frame(dom=cddids),cddid_all,by='dom',all.x=T)
}

#' @rdname Rgator-deprecated
#' @export
getDomainSets <- function( inputsites, sitecol, domaindata,  max_dom_proportion=0.81, stem_distance=100  ) {
  .Deprecated('calculateDomainSets',package='Rgator')
  calculateDomainSets(inputsites, sitecol, domaindata,  max_dom_proportion, stem_distance)
}

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
    # This is a type II transmembrane or it is secreted
    if (( (filtered$dom[1] == "SIGNALP") | grepl("tmhmm-TMhelix",filtered$dom[1]) ) & (filtered$uniprot[1] %in% typeii )) {
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