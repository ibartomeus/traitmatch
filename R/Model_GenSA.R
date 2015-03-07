############################################
#
# Bayesian model of trophic interactions
# based on body size, using observation of interactions 
# only (presence)
# 
# Dominique Gravel
# June 25th, 2014
#
############################################


######################
# The model
model = function(pars,data) {

	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]
	meanMprey = mean(unique(MPrey))
	sdMprey = sd(unique(MPrey))
	MPred = data$MPred
	MPrey = data$MPrey

	# Optimum and range
	o = a0 + a1*MPred 
	r = b0 + b1*MPred
		
	# Compute the conditional
	pLM = exp(-(o-MPrey)^2/2/r^2)

	# Compute the marginal
	pM = dnorm(x=MPrey,mean=meanMprey,sd=sdMprey)
	
	#  Integrate the denominator
	pL = r/(r^2+sdMprey^2)^0.5*exp(-(o-meanMprey)^2/2/(r^2+sdMprey^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################
