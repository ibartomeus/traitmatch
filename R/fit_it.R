#' @name fit_it
#' 
#' @title Fit a model using SA.
#' 
#' @description Fit the paramaters of a given model to data using SA. Can be slow for large datasets. 
#' 
#' @param model One of the three models implemented integrated_model, niche_model or neutral_model.
#' @param Tlevel1 Vector of trait values of the first interaction partner.
#' @param Tlevel2 Vector of trait values of the second interaction partner.
#' @param mean_Tlevel1 Mean of trait values of the first interaction partner. Can be weighted or not,
#'  and can use independent information on the trait distribution to be calculated.
#' @param sd_Tlevel1 Standard deviation of trait values of the first interaction partner. Can be weighted or not,
#'  and can use independent information on the trait distribution to be calculated.
#' @param pars a vector of the form c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) with the initial parameters.
#' @param pars_lo a vector of the form c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) with the lower limits that the 
#' parameters can reach. You can use a priosy information to constrain the posible values of the parameters 
#' (e.g. the slope, a1, has to be positive)
#' @param pars_hi a vector of the form c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) with the higher limits that the 
#' parameters can reach. You can use a priosy information to constrain the posible values of the parameters 
#' 
#' @return a vector with the estimated parameters of the form c(a0 = x, a1 = y, b0 = z, b1 = w)
#'
#' @examples 
#' See readme file here: https://github.com/ibartomeus/trait_match
#' 
#' @author
#' Dominique Gravel and Ignasi Bartomeus
#'
#' @export
fit_it <- function(model, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1,
                   pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                   par_lo = c(a0 = -5000, a1 = -5000, b0 = -5000, b1 = -5000),
                   par_hi = c(a0 = 5000, a1 = 5000, b0 = 5000, b1 = 5000)){
  require(GenSA)
  # Initial values 
  #pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) #default
  #lm_M = lm(Tlevel1 ~ Tlevel2)
  #pars = c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0) #another option
  # Maximum likelihood estimation
  estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                     control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), 
                     Tlevel1 = Tlevel1, Tlevel2 = Tlevel2, 
                     mean_Tlevel1 = mean_Tlevel1, sd_Tlevel1 = sd_Tlevel1)
  estim.pars$par  
}


