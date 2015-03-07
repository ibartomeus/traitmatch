######################
# Source packages
library(GenSA)

# Load the model
source('R/Model_GenSA.R')

# Load data
original_data = read.csv2("data/size_Barnes2008.csv",dec=".")
head(original_data)
MPred = log10(original_data$standardised_predator_length)
MPrey = log10(original_data$si_prey_length)

# Initial values 
lm_M = lm(MPrey~MPred)
pars = c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0)

# Boundaries
par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)

# Define the dat
data = data.frame(MPrey = MPrey, MPred = MPred)

# Maximum likelihood estimation
estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                   control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), data = data)

write.table(estim.pars$par,file="R/pars.txt")

