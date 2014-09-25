#'  Load a Reactome pathway network from an URL or local file
#'
#'  @param  pathwayURL URL to retrieve pathway from (http:// or file://)
#'  @return Network object for the given pathway
#'  @export
#'  @examples
#'  # Load from a local file
#'  getReactomePathway("file://mypathway.biopax")
#'  # Load from ReactomeDB
#'  url <- 'http://www.reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/535734'
#'  network <- getReactomePathway(url)
# @importFrom rBiopaxParser readBiopax
getReactomePathway <- function(pathwayURL) {
  biopax <- loadReactomePathway(pathwayURL)
  pathway2network(biopax, rBiopaxParser::listInstances(biopax, class="pathway")$id[1])
}

loadReactomePathway <- function(pathwayURL) {
  td = tempdir()
  tf = tempfile(tmpdir=td,fileext='.biopax')
  download.file(pathwayURL,tf)
  rBiopaxParser::readBiopax(tf)
}

#' Plot a Reactome network from a given pathway URL
#'
#'  @param  pathwayURL URL to retrieve pathway from (http:// or file://)
#'  @return List containing two elements, plot: Containing the plot, and coords: Containing the co-ordinates of each protein
#'  @export
generateReactomeNetwork <- function(pathwayURL) {
  ggnet(getReactomePathway(pathwayURL))
}


striphash <- function(str)
{
  gsub("#","",str)
}

#  Split up interactions in a biopax pathway so that they are their own smaller
#  sub-networks. We split a single interaction into the left and right parts
#  with a new interaction entity that sits in between the left and the right
# @importFrom rBiopaxParser getReferencedIDs listInstances
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

# Modified from pathway2regulatoryGraph in rBiopaxReader
#
# @importFrom rBiopaxParser listPathwayComponents listComplexComponents selectInstances getSubClasses getInstanceClass getXrefAnnotations
# @importFrom data.table setkeyv rbindlist
# @importFrom network network set.edge.attribute set.vertex.attribute network.vertex.names
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

#  Draw a wedge at a given point
#  @param  idx     Index of segment around the circle
#  @param  total   Total number of segments aroudn the circle
#  @param  start_radius  Inner (starting) radius
#  @param  width   Width of segment
#  @param  c_x     X Center of circle that segment is placed around
#  @param  c_y     Y Center of circle that segment is placed around
#  @param  values  data.frame with two columns minimum : datasource and value
#  @param  scales  Color scales used to color the segments
#  @seealso \code{\link{makeScale}}
make_wedges.protein <- function(idx,total,start_radius,width,c_x,c_y,values,zeroes) {
  # Arguments:  Index around circle
  #             Total segments
  #             Inner (starting radius)
  #             Width of wedge
  #             Center of wedge
  #             Center of wedge
  #             Expression values frame
  values$idx <- (0:(nrow(values)-1))
  # Rnaseq : blue -> orange
  # ms : green -> purple
  apply(values,1,function(exp_row) {
    inner_range <- seq(0,2*pi/total,length.out=30) + idx*2*pi/total
    outer_range <- rev(inner_range)
    radius <- start_radius + (0.5 + as.numeric(exp_row['idx'])) * width
    inner_line <- sapply(inner_range, function(x) { c(x=c_x+(radius - 0.5*width)*cos(x) , y=c_y+(radius - 0.5*width)*sin(x)) })
    outer_line <- sapply(outer_range, function(x) { c(x=c_x+(radius + 0.5*width)*cos(x) , y=c_y+(radius + 0.5*width)*sin(x)) })
    value <- as.numeric(exp_row['value'])
    scaledvalue <- as.numeric(exp_row['scaledvalue'])
    if (unique(as.numeric(exp_row['value']) != 0)) {
      geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',]),value=value,fillval=scaledvalue),aes(x=x,y=y,fill=fillval),color='white')
    } else {
      geom_polygon(data=data.frame(x=c(inner_line['x',],outer_line['x',]),y=c(inner_line['y',],outer_line['y',]),value=value,fillval=scaledvalue),aes(x=x,y=y),fill=zeroes[[ exp_row['datasource'] ]],color='gray')
    }
  })
}

# @importFrom scales rescale
makeScale <- function(name,in.cols,limits,fn) {
  #http://stackoverflow.com/questions/18700938/ggplot2-positive-and-negative-values-different-color-gradient
  vals <- c(min(limits),0-.Machine$double.eps,0,0+.Machine$double.eps,max(limits))
  breaks <- min(limits):max(limits)
  cols <- c(in.cols$negative[1],in.cols$negative[2],in.cols$middle,in.cols$positive[1],in.cols$positive[2])
  #name=name,
  return (list(palette=fn(colours=cols,breaks=breaks,limits=limits,values=vals)$palette, legend=fn(colours=cols,limits=limits,breaks=breaks,values=scales::rescale(vals)  )))
}

makeScaleMultiple <- function(...) {
  all.cols <- list(...)
  if (length(all.cols) < 1) {
    return (list(scale=ggplot2::scale_fill_identity(),rescaler=list(),zeroes=c() ))
  }
  order_mag <- 0
  vals <- c()
  cols <- c()
  label_scales <- list()
  for (coltype in names(all.cols)) {
    in.cols <- all.cols[[coltype]]
    limits <- in.cols$limits
    vals <- c(vals, order_mag + min(limits),order_mag-.Machine$double.eps,order_mag,order_mag+.Machine$double.eps,order_mag + max(limits) )
    cols <- c(cols, in.cols$negative[1],in.cols$negative[2],in.cols$middle,in.cols$positive[1],in.cols$positive[2])
    label_scales[[coltype]] <- order_mag
    order_mag <- order_mag + 20
  }
  scaledvals <- scales::rescale(vals)
  for (i in seq(3,length(vals),by=5)) {
    scaledvals[i-1] <- scaledvals[i] - .Machine$double.eps
    scaledvals[i+1] <- scaledvals[i] + .Machine$double.eps
  }
  list( scale=ggplot2::scale_fill_gradientn(guide="legend", colours=cols,breaks=seq(min(vals),max(vals),by=0.5),values=scaledvals,limits=c(min(vals),max(vals))) ,
        rescaler = label_scales,
        zeroes = sapply(all.cols, function(x) { x$middle },USE.NAMES=T,simplify=F) )
}

require(proto)
GeomTextBackground <- proto(ggplot2:::GeomText, {
  objname <- "backgroundtext"
  draw <- function(., data, scales, coordinates, ..., parse = FALSE,
                   expand = 1.2, bgcol = "grey50", bgfill = NA, bgalpha = 1) {
    lab <- data$label
    if (parse) {
      lab <- parse(text = lab)
    }
    with(ggplot2:::coord_transform(coordinates, data, scales), {
      sizes <- plyr::llply(1:nrow(data),
        function(i) with(data[i, ], {
          grobs <- grid::textGrob(lab[i], default.units="native", rot=angle, gp=grid::gpar(fontsize=size * .pt))
          list(w = grid::grobWidth(grobs), h = grid::grobHeight(grobs))
        }))
      widths = do.call(grid::unit.c, lapply(sizes, "[[", "w")) * expand
      heights = do.call(grid::unit.c, lapply(sizes, "[[", "h")) * expand
      rectgrobs <- as.vector(lapply(1:length(x),function(idx) { 
        grid::roundrectGrob(x[idx], y[idx],
                     width = widths[idx],
                     height = heights[idx],
                     r = heights*0.25,
                     gp = grid::gpar(col = scales::alpha(bgcol, bgalpha), fill = scales::alpha(bgfill, bgalpha)))
      }))
      grid::gList( do.call(grid::gList,rectgrobs), .super$draw(., data, scales, coordinates, ..., parse))
    })
  }
})

#' Rendering of text labels with a roundrect background on them
#'
#' Works the same as geom_text, but you can set extra parameters on the background
#' @param expand  Ratio to expand the text label background by (defaults to 1.2)
#' @param bgcol   Colour to use to draw the background roundrect
#' @param bgfill  Fill colour to use for the background roundrect
#' @param bgalpha Alpha opacity for the background roundrect
#' @export
geom_text_background <- function (...) {
  GeomTextBackground$new(...)
}

#testing_expression <- data.frame(uniprot=c(rep('Q9UHC9',2),rep('O95477',2)),value=c(3,5,-3,0),datasource=c('rnaseq','ms','rnaseq','ms'))

#' @rdname Rgator-deprecated
#' @export
overlayExpression.protein <- function(organism,plot,expression,uniprot.symbols,node.color,...) {
  networkPlot.annotate(organism,plot,expression,uniprot.symbols,node.color,...)
}

#' Annotate a protein based network plot with data from a table of annotations
#'
#' Given a network plot obtained from String \code{\link{generateStringNetwork}} or Reactome \code{\link{generateReactomeNetwork}}
#' overlay expression information as a set of halos.
#' @param organism  NCBI taxonomy id for the organism that we are working with
#' @param plot      Result from either \code{\link{generateStringNetwork}} or \code{\link{generateReactomeNetwork}}
#' @param annotations  Data table with three columns - \code{uniprot},\code{value},\code{datasource}
#' @param uniprot.symbols Toggle if the uniprot symbols should be overlaid on the graph
#' @param node.color      Datasource string value to use for setting the colour of the node
#' @param ...             Colour ranges to use for the scales for each of the given datasources
#' @export
networkPlot.annotate <- function(organism=9606,plot,annotations=data.frame(uniprot=NA,value=NA,datasource=NA),uniprot.symbols=F,node.color=NA,...) {
  coords <- plot$coords
  newplot <- plot$plot

  scale_info <- makeScaleMultiple(...)

  annotations$value <- as.numeric(annotations$value)
  annotations <- subset(annotations,!is.na(value))
  if (nrow(annotations) > 0) {
    annotations$idx <- row.names(annotations)
    annotations$size <- 1
    annotations <- subset(annotations, uniprot %in% coords$uniprot)
    annotations$scaledvalue <- apply(annotations[,c('datasource','value')],1,function(row) { as.numeric(row['value']) + scale_info$rescaler[[ row['datasource'] ]]  })
  }

#  scale_limits <- list(...)
#  scales <- list(scale1=makeScale( names(scale_limits)[1],list(negative=c('#FC7F23','#FCBE75'),middle='#ffffff',positive=c('#A7CEE2','#2579B2')), scale_limits[[1]], ggplot2::scale_color_gradientn ),
#                 scale2=makeScale( names(scale_limits)[2],list(negative=c('#399F34','#B3DE8D'),middle='#99ff99',positive=c('#CAB2D5','#6A4098')), scale_limits[[2]], ggplot2::scale_fill_gradientn )
#                 )
  #names(scales) <- names(scale_limits)
  if ( ! is.na(node.color) ) {
    node.frame <- merge(coords,subset(annotations, datasource %in% node.color),by='uniprot')
    newplot <- newplot + geom_point(data=node.frame,aes(x=X1,y=X2),size=6,color=scale_info$scale$palette(scales::rescale(node.frame$scaledvalue,from=scale_info$scale$limits)))
  }
  grobs <- apply(subset(coords,uniprot != ''),1,function(rowdata) {
    all_uniprots <- gsub("-.*","",unlist(strsplit(rowdata['uniprot']," ")))
    return_data <- c()
    for(radius in 1:length(all_uniprots)) {
      uprot <- all_uniprots[radius]
      local_expr <- subset(annotations,uniprot == uprot & ! datasource %in% node.color )
      local_expr <- local_expr[with(local_expr, order(datasource,decreasing=T)), ]
      local_expr$radius <- rep(radius,nrow(local_expr))
      if (nrow(local_expr) > 0) {
        return_data <- c(return_data, make_wedges.protein(radius,length(all_uniprots), 1, 1.5, as.numeric(rowdata['X1']), as.numeric(rowdata['X2']),local_expr,scale_info$zeroes))
      }
    }
    return (return_data)
  })
  for(grob in grobs) {
    newplot <- newplot + grob
  }
  if (uniprot.symbols) {
#    newplot <- newplot + geom_text_background(data=coords,aes(x=X1,y=X2-1),bgfill='#000000',bgalpha=0.7,expand=2,size=4,color='white',label=sapply ( strsplit(plot$coords$uniprot,' '), function(prots) { if (length(prots) > 0) { return(paste(getGeneNames(9606,prots)$symbol,collapse=' ')) } else { return("") } }  ))
    newplot <- newplot + geom_text(data=coords,aes(x=X1,y=X2-1),size=4,color='black',label=sapply ( strsplit(plot$coords$uniprot,' '), function(prots) { if (length(prots) > 0) { return(paste(getGeneNames(9606,prots)$symbol,collapse=' ')) } else { return("") } }  ))
  }
  newplot <- newplot + scale_info$scale + theme(legend.key.size = unit(0.25, "cm")) #+ scales[[1]]$legend + scales[[2]]$legend
  return(newplot)
}
