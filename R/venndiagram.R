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
  package_attached = TRUE
  browser()
  if (! 'package:VennDiagram' %in% search()) {
    require('VennDiagram')
    package_attached = FALSE
  }
  venn.plot <- generateVennDiagramGrob(data)
  if ( ! package_attached ) {
    detach('package:VennDiagram')
  }
  ggplot(data.frame(),aes(x=1,y=1)) +
    geom_blank() + ggplot2::theme_bw() + ggplot2::annotation_custom(grob = venn.plot)  +
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
  return (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(getNamespace('VennDiagram')$venn.diagram(data,filename=NULL)))), cl=c("arrange", "ggplot")))
}

