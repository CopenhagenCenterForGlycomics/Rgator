# @importFrom GO.db GOBPPARENTS GOCCPARENTS GOMFPARENTS
# @importFrom AnnotationDbi toTable
# @importFrom network set.edge.attribute set.vertex.attribute network.vertex.names

#' Create a GO network for a given set of GO identifiers
#'
#'  @param  goids     GO ids to use as the basis to build the network
#'  @param  ancestors Number of ancestors to add to the graph to establish more remote connections
#'  @return Network object with the GO network
#'  @export
getGONetwork <- function(goids,ancestors=4) {
  goNetwork(goids,ancestors)
}

goNetwork <- function(goids,ancestors=4) {
  BP <- subset(AnnotationDbi::toTable(GO.db::GOBPPARENTS), RelationshipType %in% c('is_a','part_of'))
  CC <- subset(AnnotationDbi::toTable(GO.db::GOCCPARENTS), RelationshipType %in% c('is_a','part_of'))
  MF <- subset(AnnotationDbi::toTable(GO.db::GOMFPARENTS), RelationshipType %in% c('is_a','part_of'))
  go.vertices <- rbind(BP,CC,MF)
  names(go.vertices) <- c('left','right','relationship')
  if (!is.null(goids)) {
    goids <- as.character( unique(c(goids, unlist(sapply( goids, function(x) { tail( getGOParents(x,ontology='BP'), ancestors ) } )))) )
    go.vertices <- subset(go.vertices,left %in% goids & right %in% goids)
  }
  net <- network::network( subset(as.data.frame(go.vertices),select=c(1,2)), directed=F )
  network::set.edge.attribute(net,attrname="type",value=rep(NA,nrow(go.vertices)))
  network::set.vertex.attribute(net,attrname='class',value=as.character(sapply(network::network.vertex.names(net),function(id) { "Protein" })))
  network::set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network::network.vertex.names(net),function(id) { id })))
  net
}

#' Plot a GO network for a given set of GO identifiers
#'
#'  @param  goids     GO ids to use as the basis to build the network
#'  @param  ancestors Number of ancestors to add to the graph to establish more remote connections
#'  @return List containing two elements, plot: Containing the plot, and coords: Containing the co-ordinates of each protein
#'  @export
generateGONetwork <- function(goids,ancestors=4) {
  ggnet(getGoNetwork(goids,ancestors))
}


make_wedges.go <- function(idx,total,start_radius,width,c_x,c_y,values,scales) {
  # Arguments:  Index around circle
  #             Total segments
  #             Inner (starting radius)
  #             Width of wedge
  #             Center of wedge
  #             Center of wedge
  #             Expression values list

  values$idx <- (0:(nrow(values)-1))
  # Rnaseq : blue -> orange
  # ms : green -> purple

  apply(values,1,function(exp_row) {
    inner_range <- seq(0,2*pi/total,length.out=30) + idx*2*pi/total
    outer_range <- rev(inner_range)
    radius <- start_radius + (as.numeric(exp_row['idx']) + 0.5) * width
    inner_line <- sapply(inner_range, function(x) { c(x=c_x+(radius - 0.5*width)*cos(x) , y=c_y+(radius - 0.5*width)*sin(x)) })
    outer_line <- sapply(outer_range, function(x) { c(x=c_x+(radius + 0.5*width)*cos(x) , y=c_y+(radius + 0.5*width)*sin(x)) })
    value <- exp_row['value']
    if (value %in% names(scales)) {
      color <- scales[[ value ]]
    } else {
      color <- c("white","black")
    }
    geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',]),value=c(value,value)),aes(x=x,y=y,fill=value,color=value))
  })
}

# testing_go <- data.frame(GOBPID=c('GO:0048285','GO:0000278','GO:0051301','GO:0051301','GO:0006260','GO:0000722'),name=c('set1','set1','set1','set2','set2','set2'))

#> head(goframe)
#      GOBPID       Pvalue OddsRatio  ExpCount Count Size                                   Term         name
#1 GO:0048285 8.732266e-24  4.103702 26.618175    88  534                      organelle fission rnaseq.t1.up
#2 GO:0000278 7.415948e-22  6.047615 11.588136    54  251                     mitotic cell cycle rnaseq.t1.up
#3 GO:0051301 3.011962e-20  3.554836 29.941380    88  609                          cell division rnaseq.t1.up
#4 GO:0006260 1.487761e-18  6.249035  9.199798    44  186                        DNA replication rnaseq.t1.up
#5 GO:0000722 1.256583e-16 36.767949  1.296016    17   26 telomere maintenance via recombination rnaseq.t1.up
#6 GO:0044770 1.281635e-16  4.647469 13.507184    51  275            cell cycle phase transition rnaseq.t1.up

#' @rdname Rgator-deprecated
#' @export
overlayGOValues <- function(plot,goframe,label.terms=F,...) {
  GOnetworkPlot.annotate(plot,goframe,label.terms,...)
}


#' Annotate a GO based network plot with data from a table of GO term datasource mappings
#'
#' Given a network plot obtained from  a GO network ( \code{\link{generateGONetwork}})
#' overlay annotation information as a set of halos.
#' @param plot      Result from \code{\link{generateGONetwork}}
#' @param annotations     Data table with two columns - \code{GOID},\code{datasource}. Can also use GOBPID, GOCCID and GOMFID as GO identifier column
#' @param label.terms     Toggle if the GO terms should be overlaid on the graph
#' @param ...             Colours to use for each of the given datasources
#' @export
GOnetworkPlot.annotate <- function(plot,annotations=data.frame(GOID=c(),datasource=c()),label.terms=F,...) {
  if ('GOBPID' %in% names(annotations) & ! 'GOID' %in% names(annotations)) {
    annotations$GOID <- annotations$GOBPID
  }
  if ('GOMFID' %in% names(annotations) & ! 'GOID' %in% names(annotations)) {
    annotations$GOID <- annotations$GOMFID
  }
  if ('GOCCID' %in% names(annotations) & ! 'GOID' %in% names(annotations)) {
    annotations$GOID <- annotations$GOCCID
  }

  getBiocLiteLib('GO.db')
  newplot <- plot$plot
  coords <- plot$coords
  scales <- list(...)
  if (length(scales) < 1) {
    scales <-  list(
      c("red","gray"),
      c("lightsalmon","gray"),
      c("blue","gray"),
      c("lightskyblue","gray"),
      c("forestgreen","gray"),
      c("lightgreen","gray"),
      c("white","red"),
      c("white","blue"),
      c("white","forestgreen")
    )
    minlength <- length(scales) - length( unique(annotations$datasource) )
    if (minlength < 0) {
      minlength <- 0
    }
    names(scales) <- c( head(unique(annotations$datasource),length(scales)), rep(NA, minlength) )
    scales <- scales[ which(! is.na(names(scales))) ]
  }

  grobs <- apply(subset(coords,uniprot != ''),1,function(rowdata) {
    goid <- rowdata['uniprot']
    sources <- unique(subset(annotations,GOID==goid)$datasource)
    return_data <- c()
    if (length(sources) < 1) {
      return (return_data)
    }
    for(radius in 1:length(sources)) {
        return_data <- c(return_data, make_wedges.go(radius,length(sources), 0.5,1, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),data.frame(value=sources[radius]),scales))
    }
    return (return_data)
  })
  for(grob in grobs) {
    newplot <- newplot + grob
  }
  newplot <- newplot + scale_color_manual(values=sapply(scales,function(x) { x[2] },USE.NAMES=T)) + scale_fill_manual(values=sapply(scales,function(x) { x[1] },USE.NAMES=T))
  if (label.terms) {
    mapping <- as.data.frame(AnnotationDbi::toTable(GO.db::GOTERM[coords$uniprot]))
    coords$terms <- plyr::mapvalues(coords$uniprot,from=mapping$go_id,to=mapping$Term,warn_missing=F)
    newplot <- newplot + geom_text(data=subset(coords,uniprot %in% annotations$GOID),aes(x=X1,y=X2,label=terms),size=2,hjust=0 )
  }
  return(newplot)
}