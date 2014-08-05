# @importFrom grid gTree
# @importFrom grid gList
# @importFrom grid grid.grabExpr
# @importFrom grid grid.draw
# @importFrom VennDiagram venn.diagram
#' @export
generateVennDiagram <- function(data=list(),title="Venn Diagram") {
  require(grid)
  require(gridExtra)
  return (grid::gTree(children = grid::gList(grid::grid.grabExpr(grid::grid.draw(VennDiagram::venn.diagram(data,filename=NULL,main=title)))), cl=c("arrange", "ggplot")))
}
