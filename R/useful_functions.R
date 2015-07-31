#' @name Useful functions
#' 
#' @title weighted mean and sd
#' 
#' @description just weighted mean and sd.  
#' 
#' @param x .
#' @param w
#' 
#' @return x
#'
#' @examples x
#' 
#' @export
weighted.mean <- function(x, w) { 
  sum.w <- sum(w) 
  sum(x * w) / sum(w) 
} 

weighted.sd <- function(x, w) { 
  sum.w <- sum(w) 
  sum.w2 <- sum(w^2) 
  mean.w <- sum(x * w) / sum(w) 
  ((sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2))^0.5
} 

