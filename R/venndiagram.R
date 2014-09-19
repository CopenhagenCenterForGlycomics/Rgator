# @importFrom grid gTree
# @importFrom grid gList
# @importFrom grid grid.grabExpr
# @importFrom grid grid.draw
# @importFrom VennDiagram venn.diagram
#' Generate Venn diagrams
#' @param	...		Vectors to draw a venn diagram for
#' @param	title	Title for the venn diagram
#' @return	Venn diagram plot
#' @examples
#' vennDiagram(title="Venn diagram of numbers",a=c(1,2,3),b=c(2,3,4),c=c(1,3))
#' @export
vennDiagram <- function(title="Venn diagram",...) {
  require(grid)
  require(gridExtra)
  return (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(VennDiagram::venn.diagram(list(...),filename=NULL,main=title)))), cl=c("arrange", "ggplot")))
}

generateVennDiagram <- function(data=list(),title="Venn Diagram") {
  require(grid)
  require(gridExtra)
  return (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(VennDiagram::venn.diagram(data,filename=NULL,main=title)))), cl=c("arrange", "ggplot")))
}

