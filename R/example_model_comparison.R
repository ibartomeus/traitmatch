############################################
#
# Example of model comparison using the 
# bayesian approach to fit niche models,
# taking into account abundance and the 
# distribution of traits
#
############################################

# Settings
S = 50
ID = c(1:S)

# True parameters
a0 = 1.5
a1 = 1
b0 = 1
b1 = 0.0

true_pars = c(a0 = a0, a1 = a1, b0 = b0, b1 = b1)

# Random traits
Tlevel1 = rnorm(S,1,1)
Tlevel2 = rnorm(S,0,1)

# Abundance (related to trait 1, 
# so that the distribution is not uniform)
w = ceiling(rlnorm(S,Tlevel1))

# Weight the distribution of trait 1
wTlevel1 = rep(Tlevel1,times = w)
wID = rep(ID,times = w)
eTraits = expand.grid(wTlevel1,Tlevel2)
eID = expand.grid(wID,ID)
ePairs = paste(eID[,1],eID[,2])

# Moments of the distribution of trait 1
mean_Tlevel1 = mean(Tlevel1)
mean_Tlevel1_w = weighted.mean(Tlevel1,w)

sd_Tlevel1 = sd(Tlevel1)
sd_Tlevel1_w = weighted.sd(Tlevel1,w)

# Optimum and range
o = a0 + a1*eTraits[,2] 
r = b0 + b1*eTraits[,2]

# Compute interaction probability between pairs of species
pL = exp(-(o-eTraits[,1])^2/2/r^2)

# Determine interactions
L = numeric(length(pL))
L[runif(length(pL)) < pL] = 1

# Subset the data to keep only observed interactions
# Individual-level dataset
indTlevel1 = eTraits[L!=0,1]
indTlevel2 = eTraits[L!=0,2]

# Species-level dataset
# Number of interactions per species pair
agL = numeric(S^2)
agTraits = matrix(nr = S^2, nc = 2) 

n = 1
for(i in 1:S)
  for(j in 1:S) {
      if(i!=j) agL[n] = sum(L[eID[,1]== i & eID[,2]== j])
      agTraits[n,1] = Tlevel1[i]
      agTraits[n,2] = Tlevel2[j]
      n = n+1
  }

sTlevel1 = agTraits[agL!=0,1]
sTlevel2 = agTraits[agL!=0,2]


##############################
##############################
# Model comparison
 
# Load the model
source('R/Model_GenSA.R')  

######
# Fit models
# Integrated model
res_integrated = fit_it(integrated_model, sTlevel1, sTlevel2, mean_Tlevel1, sd_Tlevel1)

# Integrated model with weights
res_integrated_w = fit_it(integrated_model, sTlevel1, sTlevel2, mean_Tlevel1_w, sd_Tlevel1_w)

# Niche model
res_niche = fit_it(niche_model, sTlevel1, sTlevel2, mean_Tlevel1, sd_Tlevel1)

######
# Compare likelihoods
# Integrated model
-integrated_model(res_integrated, sTlevel1, sTlevel2, mean_Tlevel1, sd_Tlevel1)

# Integrated model with weights
-integrated_model(res_integrated_w, sTlevel1, sTlevel2, mean_Tlevel1_w, sd_Tlevel1_w)

# Niche model
-niche_model(res_niche, sTlevel1, sTlevel2, mean_Tlevel1, sd_Tlevel1)

# Neutral model 
-neutral_model(sTlevel1, sTlevel2, mean_Tlevel1, sd_Tlevel1)

# Neutral model with weights
-neutral_model(sTlevel1, sTlevel2, mean_Tlevel1_w, sd_Tlevel1_w)









