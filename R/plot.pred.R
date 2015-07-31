#' @name plot.pred
#' 
#' @title Plot model output
#' 
#' @description .  
#' 
#' @param pars
#' @param Tlevel1
#' @param Tlevel2
#' @param xlab
#' @param ylab
#' 
#' @return x
#'
#' @examples x
#' 
#' @author
#' Dominique Gravel and Ignasi Bartomeus
#' 
#' @references
#' Williams et al 2010?
#' 
#' @export
plot.pred <- function(pars,Tlevel1, Tlevel2, xlab = "Trait level 2", ylab = "Trait level 1"){
  seqX = seq(min(Tlevel2),max(Tlevel2),0.01)
  seqY = seq(min(Tlevel1),max(Tlevel1),0.01)
  
  XY = expand.grid(seqX,seqY)
  
  # Optimum and range
  o = pars[1] + pars[2]*XY[,1]
  r = pars[3] + pars[4]*XY[,1]
  
  # Compute the conditional
  pLM = exp(-(o-XY[,2])^2/2/r^2)
  
  Z = matrix(pLM,nr = length(seqX), nc = length(seqY))
  
  image(seqX,seqY,Z,xlab = xlab ,ylab = ylab,
        col=heat.colors(10000),cex.axis = 1.25, cex.lab = 1.5, las = 1)
  points(Tlevel2,Tlevel1,pch = 19, cex = 0.5)  
}

