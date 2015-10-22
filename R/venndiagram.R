# @importFrom grid gTree
# @importFrom grid gList
# @importFrom grid grid.grabExpr
# @importFrom grid grid.draw
# @importFrom VennDiagram venn.diagram
#' Generate Venn diagrams
#' @param ...   Vectors to draw a venn diagram for
#' @return  ggplot object
#' @examples
#' generateVennDiagram(a=c(1,2,3),b=c(2,3,4),c=c(1,3))
#' @export
generateVennDiagram <- function(...) {
  data <- list(...)
  if (! is.atomic(data[[1]]) & length(data) == 1 & is.vector(data[[1]]) & is.atomic(data[[1]][[1]])) {
    data <- data[[1]]
  }
  categories = do.call(rbind,Map(function(category) {
    data.frame(category=category,value=data[[category]])
  },names(data)))
  ggplot() +
    ggplot2::theme_bw() + geom_venn(ggplot2::aes(category=category,value=value),data=categories)  +
    theme(  line = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.border=element_blank() )
}

generateVennDiagramGrob <- function(data=list()) {
  require(grid)
  require(gridExtra)
  package_attached = TRUE
  if (! 'package:VennDiagram' %in% search()) {
    require('VennDiagram')
    package_attached = FALSE
  }
  if ('package:futile.logger' %in% search()) {
    flog.threshold(futile.logger::FATAL+1,name='VennDiagramLogger')
  }
  plot = (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(getNamespace('VennDiagram')$venn.diagram(data,filename=NULL)))), cl=c("arrange")))
  if ( ! package_attached ) {
    detach('package:VennDiagram')
  }
  return (plot)
}

#' Calculate the position-weighted-matrix for a window column
#' @export
stat_peptide <- function(mapping = NULL, data = NULL, geom = "vennDiagram",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T, ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = PeptideStat,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      ...
    )
  )
}


#' msdata to peptide stat
#' @export
PeptideStat <- ggplot2::ggproto("PeptideStat", ggplot2::Stat,
                        required_aes = c('class','peptide.key','peptide','site'),
                        default_aes = ggplot2::aes(category=category,value=value),
                        compute_panel = function(data, scales) {
                          specific.peps = data[,c('class','peptide.key','peptide','site')]
                          specific.peps = plyr::ddply(specific.peps,'peptide.key',function(peps) {
                            peps$site.key = paste(sort(peps$site),sep='-')
                            peps
                          })
                          specific.peps$site = NULL
                          specific.peps$key = paste(specific.peps$peptide,specific.peps$site.key,sep='-')
                          result = specific.peps[,c('class','key')]
                          names(result) = c('category','value')
                          result
                        }
)


#' Wrap the venn diagram into a geom
#' @export
geom_venn <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity",
                          show.legend = NA, inherit.aes = FALSE,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomVennDiagram,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      ...
    )
  )
}

#' @export
GeomVennDiagram <- ggplot2::ggproto("GeomVennDiagram", ggplot2::Geom,
                        required_aes = c('category','value'),
                        default_aes = ggplot2::aes(),
                        draw_panel = function(data, scales,coord) {
                          coords <- coord$transform(data,scales)
                          cats = unique(data$category)
                          lists = sapply(cats, function(cat) {
                            data[data$category == cat,'value']
                          },simplify=F)
                          names(lists) = cats
                          generateVennDiagramGrob(lists)
                        }
)