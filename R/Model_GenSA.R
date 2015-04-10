############################################
#
# Bayesian model of trophic interactions
# based on body size, using observation of interactions 
# only (presence)
#
############################################


######################
# Full model
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

######################
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

######################
# Niche model
niche_model = function(pars, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1) {

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
	pM = 1/(max(Tlevel1)-min(Tlevel1))
	
	#  Integrate the denominator
	pL = r/sqrt(pi)
		
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################

# Useful functions
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

######################
#fit_it <- function(model, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1){
 #
 #  require(GenSA)
 #
  # Initial values 
 # pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0)
  
  # Boundaries 
 #par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
 #par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)

  # Maximum likelihood estimation
	#GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
   # control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), 
   # Tlevel1 = Tlevel1, Tlevel2 = Tlevel2, mean_Tlevel1 = mean_Tlevel1, sd_Tlevel1 = sd_Tlevel1)$par
#}



