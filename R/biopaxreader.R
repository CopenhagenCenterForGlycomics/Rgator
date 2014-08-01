# Convert BioPax files to Edges for making a network diagram

#library(rBiopaxParser)

#lipid.pax <- readBiopax('lipid.biopax')


striphash <- function(str)
{
  gsub("#","",str)
}

interaction2components <- function (biopax, id, splitComplexes = TRUE, returnIDonly = FALSE, 
                                    biopaxlevel = 3) 
{
  if ("biopax" %in% class(biopax)) {
    biopaxlevel = biopax$biopaxlevel
  }
  else if ("data.table" %in% class(biopax)) {
  }
  else {
    stop("listInteractionComponents: parameter biopax is neither biopax object nor compatible biopax data.table")
  }
  if (biopax$biopaxlevel == 2) {
    compname = c("PARTICIPANTS", "CONTROLLER", "CONTROLLED", 
                 "LEFT", "RIGHT", "COFACTOR", "PHYSICAL-ENTITY")
    if (splitComplexes) 
      compname = c(compname, "COMPONENTS")
  }
  if (biopax$biopaxlevel == 3) {
    compname = c("participant", "controller", "controlled", 
                 "left", "right", "cofactor", "product", "template")
    if (splitComplexes) 
      compname = c(compname, "component")
  }
  id = unique(striphash(id))
  comp_list_l = (unlist(sapply(intersect(compname, c("participant","controller","left","cofactor")), function(prop) {
                                                   as.character(unique(getReferencedIDs(biopax,id, recursive = TRUE, onlyFollowProperties = c(prop))))
                                                   },simplify=T,USE.NAMES=T)))
  comp_list_r = (unlist(sapply(intersect(compname, c("controlled","right","product","template")), function(prop) {
    as.character(unique(getReferencedIDs(biopax,id, recursive = TRUE, onlyFollowProperties = c(prop))))
  },simplify=T,USE.NAMES=T)))
  names(comp_list_l) <- rep("left",length(comp_list_l))
  names(comp_list_r) <- rep("right",length(comp_list_r))
  comp_list <- c( comp_list_l, comp_list_r )
  if (is.null(comp_list)) 
    return(NULL)
  comp_list = comp_list[!is.na(comp_list) & !is.null(comp_list) & 
                          nchar(comp_list) > 0]
  if (returnIDonly) 
    return(comp_list)
    #return(unique(striphash(comp_list)))
  listInstances(biopax, id = unique(striphash(comp_list)), 
                biopaxlevel = biopaxlevel)
}

pathway2network <- function (biopax, pwid,verbose=F) 
{
  if (is.null(biopax) || !("biopax" %in% class(biopax))) 
    stop("Error: pathway2RegulatoryGraph: parameter biopax has to be of class biopax.")
  biopaxlevel = biopax$biopaxlevel
  if (is.null(pwid) || !("character" %in% class(pwid))) 
    stop("Error: pathway2RegulatoryGraph: parameter pwid has to be of class character.")
  if (!require(graph)) {
    message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!", 
                  "\n"))
    return(NULL)
  }
  pwid = unique(striphash(pwid))
  pw_component_list = listPathwayComponents(biopax, pwid, returnIDonly = T)
  if (length(pw_component_list) == 0) {
    warning("Pathway seems to have no pathway components")
    return(NULL)
  }
  pw_component_list = selectInstances(biopax, id = pw_component_list, 
                                      includeReferencedInstances = TRUE, returnCopy = TRUE)
  pw_component_list$property = tolower(pw_component_list$property)
  setkeyv(pw_component_list, cols = c("id", "class", "property"))
  pw_interactions = pw_component_list[tolower(pw_component_list$class) %chin% 
                                    tolower(getSubClasses("Interaction")), 
                                  ]
  if (length(pw_interactions$id) == 0) {
    warning("warning: pathway2RegulatoryGraph: supplied graph has no interaction pathway components. Returning NULL.")
    return(NULL)
  }
  else {
    if (verbose) {
      message(paste("Found", length(unique(pw_interactions$id)), 
                    "pathway components. Putting them together..."))
    }
  }
  edges <- rbindlist(lapply(unique(pw_interactions$id), function(i) {
    interactionbits <- interaction2components(biopax,i,returnIDonly=T)
    is_target <- length(getReferencingIDs(biopax,i,onlyFollowProperties=c("participant","controller","controlled","left","right","cofactor","product","template","component"))) > 0
    pw_edges <- NA
    left_list <- interactionbits[which(names(interactionbits) == 'left')]
    right_list <- interactionbits[which(names(interactionbits) == 'right')]
    if (is_target) {
      pw_edges <- rbind( expand.grid(left_list,i),expand.grid(i,right_list)  )
      # make nodes between left interaction bits and right interaction bits passing through self
    } else{
      # make nodes between left and right bits
      pw_edges <- expand.grid(left_list,right_list)
    }
    names(pw_edges) <- c('start','end')
    pw_edges$id <- rep(i,nrow(pw_edges))
    pw_edges$class <- rep(getInstanceClass(biopax,i), nrow(pw_edges))
    prop_type <- getInstanceProperty(lipid.pax,id=i,property='controlType')
    if (is.null(prop_type)) {
      prop_type <- NA
    }
    pw_edges$type <- rep( prop_type , nrow(pw_edges) )
    pw_edges
  }))
  net <- network(subset(edges,select=c('start','end')),directed=T)
  set.edge.attribute(net,attrname="id",value=edges$id)
  set.edge.attribute(net,attrname="class",value=edges$class)
  set.edge.attribute(net,attrname="type",value=edges$type)
  
  get_proteins <- function(id) {
    uprots <- Filter(function(x) { grepl("UniProt",x) },getXrefAnnotations(biopax,id)$annotation)
    uprots <- c(uprots,unlist(lapply(listComplexComponents(biopax,id,returnIDonly=T), get_proteins)))
    uprots <- unique(gsub("-.*","",gsub("UniProt.*:","",uprots)))
    return (uprots)
  }
  set.vertex.attribute(net,attrname='class',value=as.character(sapply(network.vertex.names(net),function(id) { getInstanceClass(biopax,id) })))
  set.vertex.attribute(net,attrname='uniprot',value=lapply(sapply(network.vertex.names(net),get_proteins ),function(x) { paste(unique(x),collapse=' ') } ))
  net
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


make_wedges <- function(idx,total,start_radius,width,c_x,c_y,values,scales) {
  # Arguments:  Index around circle
  #             Total segments
  #             Inner (starting radius)
  #             Width of wedge
  #             Center of wedge
  #             Center of wedge
  #             Expression values frame
  
  values$idx <- (1:nrow(values))*width
  # Rnaseq : blue -> orange
  # ms : green -> purple
  apply(values,1,function(exp_row) {
    inner_range <- seq(0,2*pi/total,length.out=30) + idx*2*pi/total
    outer_range <- rev(inner_range)
    radius <- start_radius + 2*as.numeric(exp_row['idx']) * width
    inner_line <- sapply(inner_range, function(x) { c(x=c_x+(radius - 0.5*width)*cos(x) , y=c_y+(radius - 0.5*width)*sin(x)) })
    outer_line <- sapply(outer_range, function(x) { c(x=c_x+(radius + 0.5*width)*cos(x) , y=c_y+(radius + 0.5*width)*sin(x)) })
    value <- as.numeric(exp_row['value'])
    if ( exp_row['datasource'] == 'ms') {
      if (value > 5) {
        value <- 5
      }
      if (value < -5) {
        value <- -5
      }
    }
    if ( exp_row['datasource'] == 'rnaseq' ) {
      if (value > 7) {
        value <- 7
      }
      if (value < -7) {
        value <- -7
      }
    }
    color <- scales[[ exp_row['datasource'] ]]$palette(value)
    if ( exp_row['datasource'] == 'ms' && abs(value) < 0.00001 ) {
      color <- '#99ff99'
    }
    if ( exp_row['datasource'] == 'rnaseq' && abs(value) < 0.00001 ) {
      #color <- '#999999'
    }
    geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',]),value=value),aes(x=x,y=y),fill=color,color='gray') #+ scales[[ exp_row['datasource'] ]]
  })
}

normalise <- function(x){(x-min(x))/(max(x)-min(x))}

makeScale <- function(name,in.cols,limits,fn) {
  vals <- c(min(limits),-0.00001,-0.00001,0,0.00001,0.00001,max(limits))
  breaks <- min(limits):max(limits)
  cols <- c(in.cols$negative[1],in.cols$negative[2],in.cols$negative[2],in.cols$middle,in.cols$positive[1],in.cols$positive[1],in.cols$positive[2])
  return (list(palette=fn(colours=cols,breaks=breaks,limits=limits,values=vals)$palette, legend=fn(name=name,colours=cols,limits=limits,breaks=breaks,values=rescale(vals)  )))
}

overlayExpression <- function(plot,in_cords,expression,uniprot.symbols=F) {
  library(scales)
  cords <- in_cords
  newplot <- plot
  expression$value <- as.numeric(expression$value)
  expression <- subset(expression,!is.na(value))
  expression$idx <- row.names(expression)
  expression$shape <- sapply(expression$datasource,function(src) { if (src == "ms") { return(17) } else { return (19) } })
  expression$color <- sapply(expression$value,function(value) { if (value > 0) { return("blue") } else { return ("red") } })
  expression$size <- normalise(abs(expression$value))  #abs(scale(expression$value,center=0,scale=0.5))
  #scales <- list(rnaseq=scale_color_gradientn(colours=c('#FC7F23','#FCBE75','#FCBE75','#999999','#A7CEE2','#A7CEE2','#2579B2'),breaks=-7:7,limits=c(-7,7),values=c(-7,-0.00001,-0.00001,0,0.00001,0.00001,7) ), ms=scale_fill_gradientn(colours=c('#399F34','#B3DE8D','#B3DE8D','#99ff99','#CAB2D5','#CAB2D5','#6A4098'),breaks=-5:5,limits=c(-5,5),values=c(-5,-0.00001,-0.00001,0,0.00001,0.00001,5) )) 
  scales <- list(rnaseq=makeScale( 'rnaseq',list(negative=c('#FC7F23','#FCBE75'),middle='#ffffff',positive=c('#A7CEE2','#2579B2')), c(-7,7), scale_color_gradientn ),
                 ms=makeScale( 'ms',list(negative=c('#399F34','#B3DE8D'),middle='#99ff99',positive=c('#CAB2D5','#6A4098')), c(-5,5), scale_fill_gradientn )
                 )
  grobs <- apply(subset(cords,uniprot != ''),1,function(rowdata) {
    all_uniprots <- gsub("-.*","",unlist(strsplit(rowdata['uniprot']," ")))
    return_data <- c()
    for(radius in 1:length(all_uniprots)) {
      uprot <- all_uniprots[radius]
      local_expr <- subset(expression,uniprot == uprot)
      local_expr <- local_expr[with(local_expr, order(datasource)), ]
      local_expr$radius <- rep(radius,nrow(local_expr))
      #local_expr$x <- as.numeric(rowdata['X1'])+(1+0.5*radius)*sin( ((radius %% 2) / 10 * pi) + (as.numeric(local_expr$idx) %% 10) * 2 / 10 * pi )
      #local_expr$y <- as.numeric(rowdata['X2'])+(1+0.5*radius)*cos( ((radius %% 2) / 10 * pi) + (as.numeric(local_expr$idx) %% 10) * 2 / 10 * pi )
      if (nrow(local_expr) > 0) {
        return_data <- c(return_data, make_wedges(radius,length(all_uniprots), 0.5,0.5, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),local_expr,scales))
      }
      #return_data <- c(return_data, ( geom_path(data=circleFun(c(as.numeric(rowdata['X1']),as.numeric(rowdata['X2'])),2*(1+0.5*radius)),aes(x,y),alpha=0.2 ) ))
      #if (nrow(local_expr) > 0) {
      #  return_data <- c(return_data, ( geom_point(aes(x=x,y=y),color=local_expr$color,shape=local_expr$shape, size=local_expr$size,data=local_expr) ))
      #}
    }
    return (return_data)
  })
  for(grob in grobs) {
    newplot <- newplot + grob
  }
  if (uniprot.symbols) {
    newplot <- newplot + geom_text(data=cords,aes(x=X1,y=X2),label=getGeneNames(9606,cords$uniprot)$symbol)
  }
  newplot <- newplot + scales$rnaseq$legend + scales$ms$legend
  return(newplot)
}

stringToNetwork <- function(uniprots) {
  library(STRINGdb)
  library(plyr)
  library(network)
  library(intergraph)
  string_db <- STRINGdb$new(species=9606)
  all_string_proteins <- string_db$get_proteins()$protein_external_id
  library(org.Hs.eg.db)
  id_mappings <- select(org.Hs.eg.db,columns=c('UNIPROT','ENSEMBLPROT'),keys=uniprots,keytype='UNIPROT')
  id_mappings$ENSEMBLPROT <- sapply ( id_mappings$ENSEMBLPROT, function(x) { paste('9606.',x,sep='') })
  names(id_mappings) <- c('uniprot','STRING_id')
  id_mappings <- unique(subset(id_mappings, STRING_id %in% all_string_proteins ))
  neighbours <- string_db$get_neighbors(id_mappings$STRING_id)
  neighbour_mappings <- select(org.Hs.eg.db,columns=c('UNIPROT','ENSEMBLPROT'),keys=gsub("9606.","",neighbours),keytype='ENSEMBLPROT')
  neighbour_mappings$ENSEMBLPROT <- sapply ( neighbour_mappings$ENSEMBLPROT, function(x) { paste('9606.',x,sep='') })
  names(neighbour_mappings) <- c('STRING_id','uniprot')
  neighbour_mappings <- subset(neighbour_mappings,select=c('uniprot','STRING_id'))
  #id_mappings <- subset( unique(rbind(id_mappings,neighbour_mappings)), uniprot %in% getUniprotIds(9606))
  
  #  id_mappings <- string_db$map( data.frame(uniprot=toupper(uniprots) ),'uniprot')
  subnetwork <- asDF(string_db$get_subnetwork(id_mappings$STRING_id))
  
  subnetwork$vertexes$name <- mapvalues(subnetwork$vertexes$name,id_mappings$STRING_id,tolower(id_mappings$uniprot),warn_missing=F)
  subnetwork$edges$V1 <- mapvalues(subnetwork$edges$V1, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  subnetwork$edges$V2 <- mapvalues(subnetwork$edges$V2, subnetwork$vertexes$intergraph_id, subnetwork$vertexes$name,warn_missing=F)
  net <- network(subset(subnetwork$edges,select=c(1,2),!is.na(V1) & !is.na(V2) ),directed=F)
  #set.edge.attribute(net,attrname="id",value=edges$id)
  set.edge.attribute(net,attrname="combined",value=subnetwork$edges$combined_score)
  set.edge.attribute(net,attrname="experimental",value=subnetwork$edges$experimental)
  set.edge.attribute(net,attrname="type",value=rep(NA,nrow(subnetwork$edges)))
  set.vertex.attribute(net,attrname='class',value=as.character(sapply(network.vertex.names(net),function(id) { "Protein" })))
  set.vertex.attribute(net,attrname='uniprot',value=as.character(sapply(network.vertex.names(net),function(id) { toupper(id) })))
  return (net)
}


