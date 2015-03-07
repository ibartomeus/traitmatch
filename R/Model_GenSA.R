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
model = function(pars,Tlevel1, Tlevel2) {

	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]
	meanTlevel1 = mean(unique(Tlevel1))
	sdTlevel1 = sd(unique(Tlevel1))

	# Optimum and range
	o = a0 + a1*Tlevel2 
	r = b0 + b1*Tlevel2
		
	# Compute the conditional
	pLM = exp(-(o-Tlevel1)^2/2/r^2)

	# Compute the marginal
	pM = dnorm(x=Tlevel1,mean=meanTlevel1,sd=sdTlevel1)
	
	#  Integrate the denominator
	pL = r/(r^2+sdTlevel1^2)^0.5*exp(-(o-meanTlevel1)^2/2/(r^2+sdTlevel1^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################
