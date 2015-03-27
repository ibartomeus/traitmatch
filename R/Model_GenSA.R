############################################
#
# Bayesian model of trophic interactions
# based on body size, using observation of interactions 
# only (presence)
# 
# Dominique Gravel
#
############################################


######################
# Full model
integrated_model = function(pars,Tlevel1, Tlevel2, mean_Trait_2, sd_Trait_2) {

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
	pM = dnorm(x=Tlevel1,mean=mean_Trait_2,sd=sd_Trait_2)
	
	#  Integrate the denominator
	pL = r/(r^2+sdTlevel1^2)^0.5*exp(-(o-meanTlevel1)^2/2/(r^2+sdTlevel1^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################
# Neutral model
neutral_model = function(pars, Tlevel1, Tlevel2, mean_Trait_2, sd_Trait_2) {
	
	# Compute the posterior probability
	pML = dnorm(x=Tlevel1,mean=mean_Trait_2,sd=sd_Trait_2)
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################
# Niche model
niche_model = function(pars,Tlevel1, Tlevel2, mean_Trait_2, sd_Trait_2) {

	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]

	# Optimum and range
	o = a0 + a1*Tlevel2 
	r = b0 + b1*Tlevel2
		
	# Compute the lposterior probability
	pML = exp(-(o-Tlevel1)^2/2/r^2)/r/sqrt(pi)
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################
