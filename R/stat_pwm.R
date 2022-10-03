#' Calculate the position-weighted-matrix for a window column
#' @export
stat_pwm <- function(mapping = NULL, data = NULL, geom = "text",
                          position = "identity",
                          show.legend = NA, inherit.aes = FALSE, backFreq=list(), drop.missing=F,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = PWMFunction,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      backFreq=backFreq,
      drop.missing=drop.missing,
      ...
    )
  )
}

#' @export
PWMFunction <- ggplot2::ggproto("PWMFunction", ggplot2::Stat,
                        default_aes = ggplot2::aes(color=..label..,size=abs(4*log2(..y..))),
                        compute_panel = function(data, scales, backFreq=list(),drop.missing=F) {
#                          message("Plotting ",length(unique(data$window))," windows")
                          pwm = calculatePWM(data,'window')
                          if (drop.missing) {
                            pwm[pwm==0]<-NA
                            missing_names = apply(pwm,2,function(col) { names(col[is.na(col)]) } )
                            missing_names[[ floor(length(missing_names)/2) + 1 ]] <- c(NA)
#                            message("Missing amino acid statistics")
#                            message(paste0(capture.output(table(unlist(missing_names))), collapse = "\n"))
#                            message(paste0(capture.output(missing_names), collapse = "\n"))
                          } else {
                            pwm[pwm==0]<- 0.001
                          }
                          bval <- plyr::laply(names(backFreq),function(x) {  pwm[x,] / unlist(backFreq[[x]]) })
                          row.names(bval)<-names(backFreq)
                          bval <- bval[names(backFreq),]
                          window_size = floor( 0.5*length(dimnames(bval)[[2]]) )
                          dimnames(bval)[[2]]<- c((-1*window_size):window_size)
                          reshape2::melt(bval,varnames=c("label","x"),value.name="y")
                        }
)