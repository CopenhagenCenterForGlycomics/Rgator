#' Superimpose a function.
stat_pwm <- function(mapping = NULL, data = NULL, geom = "text",
                          position = "identity",
                          show.legend = NA, inherit.aes = FALSE, backFreq=list(), ...) {
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
      ...
    )
  )
}

#' @export
PWMFunction <- ggproto("PWMFunction", Stat,
                        default_aes = aes(color=..label..),
                        compute_panel = function(data, scales, backFreq=list(),zero=0.001) {
                          pwm = Rgator:::calculatePWM(data,'window')
                          pwm[pwm==0]<-zero
                          bval <- plyr::laply(names(backFreq),function(x) {  pwm[x,] / backFreq[[x]] })
                          row.names(bval)<-names(backFreq)
                          hydrophobicity <- c(rep('#0000ff',6),rep('#00ff66',6),rep('#000000',8))
                          names(hydrophobicity) <-  unlist(strsplit('RKDENQSGHTAPYVMCLFIW',''))
                          chemistry <- c(rep('#00ff66',7),rep('#0000ff',3),rep('#ff0000',2),rep('#000000',8))
                          names(chemistry) <- unlist(strsplit('GSTYCQNKRHDEAVLIPWFM',''))
                          bval <- bval[names(hydrophobicity),]
                          window_size = floor( 0.5*length(dimnames(bval)[[2]]) )
                          dimnames(bval)[[2]]<- c((-1*window_size):window_size)
                          reshape2::melt(bval,varnames=c("label","x"),value.name="y")
                        }
)