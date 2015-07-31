#' @name integrated_model
#' 
#' @title Bayesian model of species interactions based on traits, using observation of interactions only (presence)
#' 
#' @description .  
#' 
#' @param pars
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
#' Dominique Gravel
#'  
#' @references
#' Williams et al 2010?
#' 
#' @export
#' 
integrated_model = function(pars, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1) {
	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]
	# Optimum and range
	o = a0 + a1*Tlevel2 
	r = b0 + b1*Tlevel2		
	# Compute the conditional
	pLM = exp(-(o-Tlevel1)^2/2/r^2)
	# Compute the marginal
	pM = dnorm(x=Tlevel1, mean = mean_Tlevel1, sd = sd_Tlevel1)
	#  Integrate the denominator
	pL = r/(r^2+sd_Tlevel1^2)^0.5*exp(-(o-mean_Tlevel1)^2/2/(r^2+sd_Tlevel1^2))	
	# Compute the posterior probability
	pML = pLM*pM/pL
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues
	return(-sum(log(pML)))		
}
# Neutral model
neutral_model = function(Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1) {
	# Compute the conditional
	pLM = 1
	# Compute the marginal
	pM = dnorm(x=Tlevel1, mean = mean_Tlevel1, sd = sd_Tlevel1)
	#  Integrate the denominator
	pL = 1
	# Compute the posterior probability
	pML = pLM*pM/pL
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues
	return(-sum(log(pML)))		
}
# Niche model
niche_model = function(pars, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1) {
	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]
	# Optimum and range
	o = a0 + a1*Tlevel2  #Here we can try polimonial functions
	r = b0 + b1*Tlevel2
	# Compute the conditional
	pLM = exp(-(o-Tlevel1)^2/2/r^2)
	# Compute the marginal
	pM = 1/(max(Tlevel1)-min(Tlevel1))
	#  Integrate the denominator
	pL = r/sqrt(pi)
	# Compute the posterior probability
	pML = pLM*pM/pL
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues
	return(-sum(log(pML)))		
}


