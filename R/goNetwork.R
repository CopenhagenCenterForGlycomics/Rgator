goNetwork <- function(orig.goids) {
  library(GO.db)
  library(network)
  BP <- subset(toTable(GOBPPARENTS), RelationshipType %in% c('is_a','part_of'))
  CC <- subset(toTable(GOCCPARENTS), RelationshipType %in% c('is_a','part_of'))
  MF <- subset(toTable(GOMFPARENTS), RelationshipType %in% c('is_a','part_of'))
  go.vertices <- rbind(BP,CC,MF)
  names(go.vertices) <- c('left','right','relationship')
  goids <- orig.goids
  if (!is.na(goids)) {
    goids <- as.character( unique(c(goids, unlist(sapply( goids, function(x) { tail( GO_parents(x,ontology='BP'), 4 ) } )))) )
    go.vertices <- subset(go.vertices,left %in% goids & right %in% goids)
  }
  net <- network( subset(as.data.frame(go.vertices),select=c(1,2)), directed=F )
  set.edge.attribute(net,attrname="type",value=rep(NA,nrow(go.vertices)))
  set.vertex.attribute(net,attrname='class',value=as.character(sapply(network.vertex.names(net),function(id) { "Protein" })))
  set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network.vertex.names(net),function(id) { id })))
  #data.frame(goid=network.vertex.names(net))
  net
}

make_wedges.go <- function(idx,total,start_radius,width,c_x,c_y,values) {
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

overlayGOValues <- function(plot,goframe,label.terms=F) {
  newplot <- plot$plot
  coords <- plot$cords
  grobs <- apply(subset(coords,uniprot != ''),1,function(rowdata) {
    goid <- rowdata['uniprot']
    sources <- unique(subset(goframe,GOBPID==goid)$name)
    return_data <- c()
    if (length(sources) < 1) {
      return (return_data)
    }
    for(radius in 1:length(sources)) {
        return_data <- c(return_data, make_wedges.go(radius,length(sources), 0,1, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),data.frame(value=sources[radius]) ))
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