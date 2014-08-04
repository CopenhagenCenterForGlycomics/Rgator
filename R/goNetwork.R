#' @importFrom Go.db GOBPPARENTS,GOCCPARENTS, GOMFPARENTS
#' @importFrom AnnotationDbi toTable
#' @importFrom network set.edge.attribute, set.vertex.attribute, network.vertex.names
goNetwork <- function(goids) {
  BP <- subset(AnnotationDbi::toTable(GO.db::GOBPPARENTS), RelationshipType %in% c('is_a','part_of'))
  CC <- subset(AnnotationDbi::toTable(GO.db::GOCCPARENTS), RelationshipType %in% c('is_a','part_of'))
  MF <- subset(AnnotationDbi::toTable(GO.db::GOMFPARENTS), RelationshipType %in% c('is_a','part_of'))
  go.vertices <- rbind(BP,CC,MF)
  names(go.vertices) <- c('left','right','relationship')
  if (!is.na(goids)) {
    goids <- as.character( unique(c(goids, unlist(sapply( goids, function(x) { tail( getGOParents(x,ontology='BP'), 4 ) } )))) )
    go.vertices <- subset(go.vertices,left %in% goids & right %in% goids)
  }
  net <- network::network( subset(as.data.frame(go.vertices),select=c(1,2)), directed=F )
  network::set.edge.attribute(net,attrname="type",value=rep(NA,nrow(go.vertices)))
  network::set.vertex.attribute(net,attrname='class',value=as.character(sapply(network::network.vertex.names(net),function(id) { "Protein" })))
  network::set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network::network.vertex.names(net),function(id) { id })))
  net
}

make_wedges.go <- function(idx,total,start_radius,width,c_x,c_y,values,scales) {
  # Arguments:  Index around circle
  #             Total segments
  #             Inner (starting radius)
  #             Width of wedge
  #             Center of wedge
  #             Center of wedge
  #             Expression values list
  
  values$idx <- (1:nrow(values))*width
  # Rnaseq : blue -> orange
  # ms : green -> purple
  if (length(scales) < 1) {
    scales <-  list(
      rnaseq.t1.up=c("red","gray"),
      diff.t1.down=c("lightsalmon","gray"),
      rnaseq.t2.up=c("blue","gray"),
      diff.t2.down=c("lightskyblue","gray"),
      rnaseq.t3.down=c("forestgreen","gray"),
      diff.t3.up=c("lightgreen","gray"),
      rnaseq.t1.down=c("white","red"),
      rnaseq.t2.down=c("white","blue"),
      rnaseq.t3.up=c("white","forestgreen")
    )
  }
  apply(values,1,function(exp_row) {
    inner_range <- seq(0,2*pi/total,length.out=30) + idx*2*pi/total
    outer_range <- rev(inner_range)
    radius <- start_radius + 2*as.numeric(exp_row['idx']) * width
    inner_line <- sapply(inner_range, function(x) { c(x=c_x+(radius - 0.5*width)*cos(x) , y=c_y+(radius - 0.5*width)*sin(x)) })
    outer_line <- sapply(outer_range, function(x) { c(x=c_x+(radius + 0.5*width)*cos(x) , y=c_y+(radius + 0.5*width)*sin(x)) })
    value <- exp_row['value']
    color <- scales[[ value ]]
    geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',])),aes(x=x,y=y),fill=color[1],color=color[2])
  })
}

#> head(goframe)
#      GOBPID       Pvalue OddsRatio  ExpCount Count Size                                   Term         name
#1 GO:0048285 8.732266e-24  4.103702 26.618175    88  534                      organelle fission rnaseq.t1.up
#2 GO:0000278 7.415948e-22  6.047615 11.588136    54  251                     mitotic cell cycle rnaseq.t1.up
#3 GO:0051301 3.011962e-20  3.554836 29.941380    88  609                          cell division rnaseq.t1.up
#4 GO:0006260 1.487761e-18  6.249035  9.199798    44  186                        DNA replication rnaseq.t1.up
#5 GO:0000722 1.256583e-16 36.767949  1.296016    17   26 telomere maintenance via recombination rnaseq.t1.up
#6 GO:0044770 1.281635e-16  4.647469 13.507184    51  275            cell cycle phase transition rnaseq.t1.up


overlayGOValues <- function(plot,goframe,label.terms=F,...) {
  newplot <- plot$plot
  coords <- plot$cords
  scales <- list(...)
  grobs <- apply(subset(coords,uniprot != ''),1,function(rowdata) {
    goid <- rowdata['uniprot']
    sources <- unique(subset(goframe,GOBPID==goid)$name)
    return_data <- c()
    if (length(sources) < 1) {
      return (return_data)
    }
    for(radius in 1:length(sources)) {
        return_data <- c(return_data, make_wedges.go(radius,length(sources), 0,1, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),data.frame(value=sources[radius]),scales))
    }
    return (return_data)
  })
  for(grob in grobs) {
    newplot <- newplot + grob
  }
  if (label.terms) {
    mapping <- as.data.frame(toTable(GOTERM[coords$uniprot]))
    coords$terms <- mapvalues(coords$uniprot,from=mapping$go_id,to=mapping$Term,warn_missing=F)
    newplot <- newplot + geom_text(data=subset(coords,uniprot %in% goframe$GOBPID),aes(x=X1,y=X2,label=terms),size=2,hjust=0 )
  }
  newplot <- newplot + scale_fill_identity()
  return(newplot)
}