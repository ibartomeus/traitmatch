######################

#fitting function----
fit_it <- function(model, Tlevel1, Tlevel2, mean_Tlevel1, sd_Tlevel1){
  require(GenSA)
  
  # Load the models
  source('R/Model_GenSA.R')
  
  # Initial values 
  pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0) #a0 and b1 as average and range of trait !!
  #lm_M = lm(Tlevel1 ~ Tlevel2)
  #pars = c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0)
  
  # Boundaries 
  #par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10) #try different constrains, Here a1 is set to be positive.
  par_lo = c(a0 = -1000, a1 = -1000, b0 = -100, b1 = -1000) #try different constrains, Here a1 is set to be positive.
  par_hi = c(a0 = 100, a1 = 5000, b0 = 100, b1 = 1000)
  #Nacho comment: are those especific to the data ranges?
  
  # Maximum likelihood estimation
  estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                     control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), 
                     Tlevel1 = Tlevel1, Tlevel2 = Tlevel2, 
                     mean_Tlevel1 = mean_Tlevel1, sd_Tlevel1 = sd_Tlevel1)
  
  estim.pars$par  
}

#ploting function
plot.pred <- function(pars,Tlevel1, Tlevel2, xlab = "Trait level 2", ylab = "Trait level 1"){
  seqX = seq(min(Tlevel2),max(Tlevel2),0.01)
  seqY = seq(min(Tlevel1),max(Tlevel1),0.01)
  
  XY = expand.grid(seqX,seqY)
  
  # Optimum and range
  o = pars[1] + pars[2]*XY[,1]
  r = pars[3] + pars[4]*XY[,1]
  
  # Compute the conditional
  pLM = exp(-(o-XY[,2])^2/2/r^2)
  
  Z = matrix(pLM,nr = length(seqX), nc = length(seqY))
  
  image(seqX,seqY,Z,xlab = xlab ,ylab = ylab,
        col=heat.colors(10000),cex.axis = 1.25, cex.lab = 1.5, las = 1)
  points(Tlevel2,Tlevel1,pch = 19, cex = 0.5)  
}


