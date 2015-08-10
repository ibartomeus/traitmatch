#' @name models
#' @aliases integrated_model
#' @aliases niche_model
#' @aliases neutral_model
#' 
#' @title Bayesian models for characterizing trait-based species interactions.
#' 
#' @description This models are described in Bartomeus at al. (2015 Functional Ecology) and are based
#' on Williams et al. 2000 models. They are intended to be used with fit_it function. The integrated model 
#' considers both neutral and niche constraints, while the neutral and niche models only consider its respective
#' components.
#' 
#' @author
#' Dominique Gravel
#'  
#' @references
#'Williams, R.J., Anandanadesan, A. & Purves, D. (2010) The probabilistic niche model reveals the niche structure and role of body size in a complex food web. PloS one, 5, e12092.
#'Williams, R.J. & Martinez, N.D. (2000) Simple rules yield complex food webs. Nature, 404, 180–183.#' 
#'
#' @rdname models 
#' @export 
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
#' @export
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
#' @export
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


