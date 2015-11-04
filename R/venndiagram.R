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
  plot = ( getNamespace('VennDiagram')$venn.diagram(data,filename=NULL) )
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
    stat = MsdataStat,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      level="peptide",
      ...
    )
  )
}

#' Calculate the position-weighted-matrix for a window column
#' @export
stat_site <- function(mapping = NULL, data = NULL, geom = "vennDiagram",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T, ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = MsdataStat,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      level="site",
      ...
    )
  )
}

#' Calculate the position-weighted-matrix for a window column
#' @export
stat_protein <- function(mapping = NULL, data = NULL, geom = "vennDiagram",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T, ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = MsdataStat,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      level="protein",
      ...
    )
  )
}


#' msdata to peptide stat
#' @export
MsdataStat <- ggplot2::ggproto("MsdataStat", ggplot2::Stat,
                        required_aes = c('class','peptide.key', 'uniprot', 'peptide', 'site'),
                        default_aes = ggplot2::aes(category=category,value=value),
                        compute_panel = function(data, scales, level=c("site","peptide","protein")) {
                          if (level == "peptide") {
                            specific.peps = data[,c('class','peptide.key','peptide','site')]
                            get_sites = function(peps) {
                              peps$site.key = paste(sort(peps$site),collapse='-')
                              peps
                            }
                            specific.peps = dplyr::do(dplyr::group_by(specific.peps,peptide.key), get_sites(.))
                            specific.peps$site = NULL
                            specific.peps$key = paste(specific.peps$peptide,specific.peps$site.key,sep='-')
                            result = specific.peps[,c('class','key')]
                          }
                          if (level == "site") {
                            specific.sites = data[,c('class','uniprot','site')]
                            specific.sites$key = paste(specific.sites$uniprot,specific.sites$site,sep='-')
                            result = specific.sites[,c('class','key')]
                          }
                          if (level == "protein") {
                            specific.prots = data[,c('class','uniprot')]
                            specific.prots$key = specific.prots$uniprot
                            result = specific.prots[,c('class','key')]
                          }
                          names(result) = c('category','value')
                          message(nrow(unique(result)))
                          unique(result)
                        }
)


#' Wrap the venn diagram into a geom
#' @export
geom_venn <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity",
                          show.legend = NA, inherit.aes = FALSE,na.rm=T,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomVennDiagram,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
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
                          cats = sort(unique(data$category))
                          lists = sapply(cats, function(cat) {
                            data[data$category == cat,'value']
                          },simplify=F)
                          names(lists) = cats
                          return(generateVennDiagramGrob(lists))
                        }
)