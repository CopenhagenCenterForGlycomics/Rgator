# Convert BioPax files to Edges for making a network diagram

#' @importFrom rBiopaxParser readBiopax
#' @export
getReactomePathway <- function(pathwayURL) {
  td = tempdir()
  tf = tempfile(tmpdir=td,fileext='.biopax')
  download.file(pathwayURL,tf)
  rBiopaxParser::readBiopax(tf)
}


striphash <- function(str)
{
  gsub("#","",str)
}

#' @importFrom rBiopaxParser getReferencedIDs, listInstances
interaction2components <- function (biopax, id, splitComplexes = TRUE, returnIDonly = FALSE,
                                    biopaxlevel = 3)
{
  if ("biopax" %in% class(biopax)) {
    biopaxlevel = biopax$biopaxlevel
  }
  if (biopax$biopaxlevel == 2) {
    compname = c("PARTICIPANTS", "CONTROLLER", "CONTROLLED",
                 "LEFT", "RIGHT", "COFACTOR", "PHYSICAL-ENTITY")
    if (splitComplexes)
      compname = c(compname, "COMPONENTS")
    leftBits <- c('PARTICIPANTS','CONTROLLER','LEFT','COFACTOR')
    rightBits <- c('CONTROLLED','RIGHT')

  }
  if (biopax$biopaxlevel == 3) {
    compname = c("participant", "controller", "controlled",
                 "left", "right", "cofactor", "product", "template")
    if (splitComplexes)
      compname = c(compname, "component")

    leftBits <- c("participant","controller","left","cofactor")
    rightBits <- c("controlled","right","product","template")
  }
  id = unique(striphash(id))
  comp_list_l = (unlist(sapply(intersect(compname, leftBits), function(prop) {
                                                   as.character(unique(rBiopaxParser::getReferencedIDs(biopax,id, recursive = TRUE, onlyFollowProperties = c(prop))))
                                                   },simplify=T,USE.NAMES=T)))
  comp_list_r = (unlist(sapply(intersect(compname, rightBits), function(prop) {
    as.character(unique(rBiopaxParser::getReferencedIDs(biopax,id, recursive = TRUE, onlyFollowProperties = c(prop))))
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
  rBiopaxParser::listInstances(biopax, id = unique(striphash(comp_list)),
                biopaxlevel = biopaxlevel)
}

#' Modified from pathway2regulatoryGraph in rBiopaxReader
#
#' @importFrom rBiopaxParser listPathwayComponents, listComplexComponents, selectInstances, getSubClasses, getInstanceClass, getXrefAnnotations
#' @importFrom data.table setkeyv, rbindlist
#' @importFrom network network, set.edge.attribute, set.vertex.attribute, network.vertex.names
pathway2network <- function (biopax, pwid,verbose=F) {
  biopaxlevel = biopax$biopaxlevel
  pwid = unique(striphash(pwid))
  pw_component_list = rBiopaxParser::listPathwayComponents(biopax, pwid, returnIDonly = T)
  if (length(pw_component_list) == 0) {
    warning("Pathway seems to have no pathway components")
    return(NULL)
  }
  pw_component_list = rBiopaxParser::selectInstances(biopax, id = pw_component_list,
                                      includeReferencedInstances = TRUE, returnCopy = TRUE)
  pw_component_list$property = tolower(pw_component_list$property)

  data.table::setkeyv(pw_component_list, cols = c("id", "class", "property"))
  pw_interactions = pw_component_list[ ! is.na( data.table::chmatch( tolower(pw_component_list$class) , tolower(rBiopaxParser::getSubClasses("Interaction"))) ),]
  if (length(pw_interactions$id) == 0) {
    warning("warning: pathway2network: supplied graph has no interaction pathway components. Returning NULL.")
    return(NULL)
  }
  else {
    if (verbose) {
      message(paste("Found", length(unique(pw_interactions$id)),
                    "pathway components. Putting them together..."))
    }
  }
  edges <- data.table::rbindlist(lapply(unique(pw_interactions$id), function(i) {
    interactionbits <- interaction2components(biopax,i,returnIDonly=T)
    is_target <- length(rBiopaxParser::getReferencingIDs(biopax,i,onlyFollowProperties=c("participant","controller","controlled","left","right","cofactor","product","template","component"))) > 0
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
    pw_edges$class <- rep(rBiopaxParser::getInstanceClass(biopax,i), nrow(pw_edges))
    prop_type <- rBiopaxParser::getInstanceProperty(biopax,id=i,property='controlType')
    if (is.null(prop_type)) {
      prop_type <- NA
    }
    pw_edges$type <- rep( prop_type , nrow(pw_edges) )
    pw_edges
  }))
  net <- network::network(subset(edges,select=c('start','end')),directed=T)
  network::set.edge.attribute(net,attrname="id",value=edges$id)
  network::set.edge.attribute(net,attrname="class",value=edges$class)
  network::set.edge.attribute(net,attrname="type",value=edges$type)

  get_proteins <- function(id) {
    uprots <- Filter(function(x) { grepl("UniProt",x) },rBiopaxParser::getXrefAnnotations(biopax,id)$annotation)
    uprots <- c(uprots,unlist(lapply(rBiopaxParser::listComplexComponents(biopax,id,returnIDonly=T), get_proteins)))
    uprots <- unique(gsub("-.*","",gsub("UniProt.*:","",uprots)))
    return (uprots)
  }
  network::set.vertex.attribute(net,attrname='class',value=as.character(sapply(network::network.vertex.names(net),function(id) { rBiopaxParser::getInstanceClass(biopax,id) })))
  network::set.vertex.attribute(net,attrname='uniprot',value=lapply(sapply(network::network.vertex.names(net),get_proteins ),function(x) { paste(unique(x),collapse=' ') } ))
  net
}

make_wedges.protein <- function(idx,total,start_radius,width,c_x,c_y,values,scales) {
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
    if ( exp_row['datasource'] %in% names(scales) ) {
      scale_limits <- scales[[ exp_row['datasource'] ]]$legend$limits
      if (value > scale_limits[2]) {
        value <- scale_limits[2]
      }
      if (value < scale_limits[1]) {
        value <- scale_limits[1]
      }
    }
    if ( abs(value) < 0.00001 ) {
      value <- 0
    }
    color <- scales[[ exp_row['datasource'] ]]$palette(value)
    geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',]),value=value),aes(x=x,y=y),fill=color,color='gray')
  })
}

#' @importFrom scales rescale
makeScale <- function(name,in.cols,limits,fn) {
  vals <- c(min(limits),-0.00001,-0.00001,0,0.00001,0.00001,max(limits))
  breaks <- min(limits):max(limits)
  cols <- c(in.cols$negative[1],in.cols$negative[2],in.cols$negative[2],in.cols$middle,in.cols$positive[1],in.cols$positive[1],in.cols$positive[2])
  return (list(palette=fn(colours=cols,breaks=breaks,limits=limits,values=vals)$palette, legend=fn(name=name,colours=cols,limits=limits,breaks=breaks,values=scales::rescale(vals)  )))
}

#testing_expression <- data.frame(uniprot=c(rep('Q9UHC9',2),rep('O95477',2)),value=c(3,5,-3,0),datasource=c('rnaseq','ms','rnaseq','ms'))

overlayExpression.protein <- function(organism=9606,plot,expression,uniprot.symbols=F,...) {
  cords <- plot$cords
  newplot <- plot$plot
  expression$value <- as.numeric(expression$value)
  expression <- subset(expression,!is.na(value))
  expression$idx <- row.names(expression)
  expression$size <- 1
  scale_limits <- list(...)
  scales <- list(scale1=makeScale( names(scale_limits)[1],list(negative=c('#FC7F23','#FCBE75'),middle='#ffffff',positive=c('#A7CEE2','#2579B2')), scale_limits[[1]], ggplot2::scale_color_gradientn ),
                 scale2=makeScale( names(scale_limits)[2],list(negative=c('#399F34','#B3DE8D'),middle='#99ff99',positive=c('#CAB2D5','#6A4098')), scale_limits[[2]], ggplot2::scale_fill_gradientn )
                 )
  names(scales) <- names(scale_limits)
  grobs <- apply(subset(cords,uniprot != ''),1,function(rowdata) {
    all_uniprots <- gsub("-.*","",unlist(strsplit(rowdata['uniprot']," ")))
    return_data <- c()
    for(radius in 1:length(all_uniprots)) {
      uprot <- all_uniprots[radius]
      local_expr <- subset(expression,uniprot == uprot)
      local_expr <- local_expr[with(local_expr, order(datasource)), ]
      local_expr$radius <- rep(radius,nrow(local_expr))
      if (nrow(local_expr) > 0) {
        return_data <- c(return_data, make_wedges.protein(radius,length(all_uniprots), 0.5,0.5, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),local_expr,scales))
      }
    }
    return (return_data)
  })
  for(grob in grobs) {
    newplot <- newplot + grob
  }
  if (uniprot.symbols) {
    newplot <- newplot + geom_text(data=cords,aes(x=X1,y=X2),label=sapply ( strsplit(reactome.plot$cords$uniprot,' '), function(prots) { if (length(prots) > 0) { return(paste(getGeneNames(9606,prots)$symbol,collapse=' ')) } else { return("") } }  ))
  }
  newplot <- newplot + scales[[1]]$legend + scales[[2]]$legend
  return(newplot)
}
