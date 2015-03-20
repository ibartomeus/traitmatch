######################

#fitting function----
fit_it <- function(Tlevel1, Tlevel2){
  require(GenSA)
  
  # Load the model
  source('R/Model_GenSA.R')
  
  # Initial values 
  lm_M = lm(Tlevel1 ~ Tlevel2)
  pars = c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0)
  
  # Boundaries 
  par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
  par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)
  #Nacho comment: are those especific to the data ranges?
  
  # Maximum likelihood estimation
  estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                     control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), 
                     Tlevel1 = Tlevel1, Tlevel2 = Tlevel2)
  
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
  points(Tlevel2,Tlevel1,pch = 19, cex = 0.1)  
}


# Load data for Pred/Prey
#I DIDN'T UPLOAD THE DATA TO GIT BECAUSE THE REPO WILL GO OPEN EVENTUALLY AND I DON'T KNOW IF THIS DATA IS PUBLIC
original_data = read.csv2("temp/size_Barnes2008.csv",dec=".")
head(original_data)
MPred = log10(original_data$standardised_predator_length)
MPrey = log10(original_data$si_prey_length)


#pred/prey----
pars_pre <- fit_it(Tlevel1 = MPred, 
                     Tlevel2 = MPrey)

plot.pred(pars = pars_pre, Tlevel1 = MPred, 
          Tlevel2 = MPrey)


# Load data for grasshoppers
grass = read.table("data/Deraison2014.txt", h = TRUE)
head(grass)
grass <- subset(grass, Herbivory > 10)

#grasshoppers--------
pars_grass_bin <- fit_it(Tlevel1 = grass$G.Incisive.strength, 
                     Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio)

plot.pred(pars = pars_grass_bin, Tlevel1 = grass$G.Incisive.strength, 
          Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio, xlab = "C:N ratio", ylab = "Incisive strength")

# Make frequentistic grasshoper data
grass = read.table("data/Deraison2014.txt", h = TRUE)
grass <- subset(grass, Herbivory > 0)
head(grass)

#NOTE: If model is sentistive to sample size, which may not be, 
  #this may be an arbitrary way to convert to frequencies.
Incisive.strength <- c()
for(i in 1:nrow(grass)){
  temp <- rep(grass$G.Incisive.strength[i], round(grass$Herbivory[i]))
  Incisive.strength <- c(Incisive.strength,temp)
}

Carbon.nitrogen <- c()
for(i in 1:nrow(grass)){
  temp <- rep(grass$P.Leaf.carbon.nitrogen.ratio[i], round(grass$Herbivory[i]))
  Carbon.nitrogen <- c(Carbon.nitrogen,temp)
}

#and fit
pars_grass_freq <- fit_it(Tlevel1 = Incisive.strength, 
                     Tlevel2 = Carbon.nitrogen)

plot.pred(pars = pars_grass_freq, Tlevel1 = grass$G.Incisive.strength, 
          Tlevel2 = grass$P.Leaf.dry.matter.content, xlab = "C:N ratio", ylab = "Incisive strength")
#comparision of bin and freq plots looks beatiful!


#Host/para-------
host = read.table("data/Tylianakis2008.txt", h = TRUE)
head(host)

plot(host$host_body_length ~ host$parasite_body_length)

pars_host_bin <- fit_it(Tlevel1 = host$host_body_length, 
                         Tlevel2 = host$parasite_body_length)

plot.pred(pars = pars_host_bin, Tlevel1 = host$host_body_length, 
          Tlevel2 = host$parasite_body_length)

#log?
pars_host_bin_log <- fit_it(Tlevel1 = log(host$host_body_length), 
                        Tlevel2 = log(host$parasite_body_length))

plot.pred(pars = pars_host_bin_log, Tlevel1 = log(host$host_body_length), 
          Tlevel2 = log(host$parasite_body_length))
#same result... not biologically meaningful.

#pollintors
pols = read.table("data/Bartomeus2008.txt", h = TRUE)
head(pols)
pols <- subset(pols, nectar_holder_depth_mm < 15)
#Dominique, note that this outyier has a tremendous influence in the model.

pars_pols <- fit_it(Tlevel1 = pols$IT_mm, 
                        Tlevel2 = pols$nectar_holder_depth_mm)

plot.pred(pars = pars_pols, Tlevel1 = pols$IT_mm, 
          Tlevel2 = pols$nectar_holder_depth_mm, xlab = "Nectar depth", ylab = "Pollinator size")




