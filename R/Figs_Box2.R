
#source("FitModel.R")

# Figure layout

quartz(width = 10, height = 3.)

par(mar = c(5,6,2,1),mfcol = c(1,3))

# Figure 1A
# Generate random traits
S = 25
T1 = runif(S,0,1)
T2 = runif(S,0,1)

# Set boundaries
lowT1 = 0.5
lowT2 = 0.4
highT1 = 0.8
highT2 = 0.8

interactions = numeric(S)+1
interactions[T1>lowT1 & T1 <highT1 & T2>lowT2 & T2<highT2] = 19

plot(T1,T2,pch = interactions, xlab = "Trait of species i", ylab = "Trait of species j",cex.axis = 1.25, cex.lab = 1.5)
abline(v = lowT1, lty = 3)
abline(v = highT1, lty = 3)
abline(h = lowT2, lty = 3)
abline(h = highT2, lty = 3)

# Figure 1B
# Subset the preys of Thunnus albacares
#target_sp = "Dicentrarchus labrax"
target_sp = "Nototheniops larseni"
target_preys = log10(original_data$si_prey_length[original_data$predator==target_sp])
hist_target = hist(target_preys,plot = FALSE, breaks = 10)
hist_target$counts = hist_target$counts/length(target_preys)
hist_preys = hist(MPrey,plot = FALSE,breaks=40)
hist_preys$counts = hist_preys$counts/length(MPrey)

plot(hist_target,col=rgb(0,0,0,1/2), xlim = range(MPrey),ylim = c(0,0.4),ylab = "Density",cex.axis = 1.25, cex.lab = 1.5,xlab = "log10 Body Size (cm)",main="")
plot(hist_preys,col=rgb(1,1,1,1/2), add = TRUE)

target_size = mean(log10(original_data$standardised_predator_length[original_data$predator == target_sp]))
pars = read.table("pars.txt")
o = pars[1,] + pars[2,]*target_size
r = pars[3,] + pars[4,]*target_size
seqM = seq(min(MPrey),max(MPrey),0.01)
pLM = exp(-(o-seqM)^2/2/r^2)
lines(seqM,0.2*pLM,lwd = 2)

# Figure 1C
seqX = seq(min(MPred),max(MPred),0.01)
seqY = seq(min(MPrey),max(MPrey),0.01)

XY = expand.grid(seqX,seqY)

# Optimum and range
o = pars[1,] + pars[2,]*XY[,1]
r = pars[3,] + pars[4,]*XY[,1]
		
# Compute the conditional
pLM = exp(-(o-XY[,2])^2/2/r^2)

Z = matrix(pLM,nr = length(seqX), nc = length(seqY))

image(seqX,seqY,Z,xlab = "log10 Predator Body Size (cm)",ylab = "log10 Prey Body Size (cm)",col=heat.colors(10000),cex.axis = 1.25, cex.lab = 1.5)
points(MPred,MPrey,pch = 19, cex = 0.1)

dev.copy2pdf(file = "fig_box2.pdf")




