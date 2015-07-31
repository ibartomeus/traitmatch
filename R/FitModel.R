#' @name fit_it
#' 
#' @title Fit a model using SA
#' 
#' @description .  
#' 
#' @param model
#' @param Tlevel1
#' @param Tlevel2
#' @param mean_Tlevel1
#' @param sd_Tlevel1
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
#' 
fit_it <- function(model, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1){
  require(GenSA)
  # Initial values 
  pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) #a0 and b1 as average and range of trait !!
  #lm_M = lm(Tlevel1 ~ Tlevel2)
  #pars = c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0)
  # Boundaries 
  #par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10) #try different constrains, Here a1 is set to be positive.
  par_lo = c(a0 = -1000, a1 = -5000, b0 = -1000, b1 = -1000) #try different constrains, Here a1 is set to be positive.
  par_hi = c(a0 = 5000, a1 = 5000, b0 = 1000, b1 = 1000)
  #Nacho comment: are those especific to the data ranges?
  # Maximum likelihood estimation
  estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                     control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), 
                     Tlevel1 = Tlevel1, Tlevel2 = Tlevel2, 
                     mean_Tlevel1 = mean_Tlevel1, sd_Tlevel1 = sd_Tlevel1)
  estim.pars$par  
}


