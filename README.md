# traitmatch: Package to predict trait matching from species interactions.

This document reprodcue the analysis done in Bartomeus et al 2015 (Functional Ecology) as a way to show how it works.

To install the package run (only once):

```
install.packages("devtools") 
require(devtools)
install_github("traitmatch", "ibartomeus")
require(traitmatch)
```

We first use data on predator and prey body size from Barnes et al. 2008. We load the data and logtransform it.

```
# Load data for Pred/Prey
fish <- read.table("data/Barnes2008.txt", header = TRUE)
#each row represents an interaction
head(fish) 
#Note that we had issues if we directly log transform the data.
#MPred = log10(fish$standardised_predator_length)
#MPrey = log10(fish$si_prey_length)
#It cause problems with the precision of the integral of the normal distribution. 
#It has very strange behaviour when the standard deviation is small (<<1). 
#We solved this by *10

MPred = log10(fish$standardised_predator_length*10)
MPrey = log10(fish$si_prey_length*10)
```

Then we fit the `integrated_model` that integrates neutral and niche constrains. We use pairwise Prey-Predatory interactions (Tlevel vectors) and we asume that the distribution of preys is well characterized by this network (we use mean and sd directly from the distribution of prey interactions). Here we constrain the parameters with a priory information. For example, the slope of the relationship has to be positive (the bigger the predator, the bigger the prey). Tuning the a priori assumptions is important to get good estimates of the parameters. With unconstrained estimates, some models, like this one, can't find an optima, and show a flat vertical slope with no warning, so always check that your results make sense before accepting them as the truth. This process is slow, specially for large datasets like this one. You can use the max.time (in seconds) to cut the process to e.g. 900 secoonds (15 minuts). The default is 30 minuts. Note that the estimates can be bad if too little time allowed. 

```
?fit_it
pars_pre <- fit_it(integrated_model, 
                   Tlevel1 = MPrey,  
                   Tlevel2 = MPred,
                   mean_Tlevel1 = mean(MPrey),
                   sd_Tlevel1 = sd(MPrey),
                   pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                   par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
                   par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
                   max.time = 900)
pars_pre 
```

We can plot the predicted model:

```
#With this model, the first number gives you the intercept and second number gives you the slope.
#plot(MPrey ~ MPred)
#abline(a = pars_pre[1], b = pars_pre[2], col = "red")
#but to get the ranges, we have a nicer function:

?plot_pred
plot_pred(pars = pars_pre, Tlevel1 = MPrey, 
          Tlevel2 = MPred, xlab = "log (Predator body size)", 
          ylab = "log (prey body size)", pch = ".")
```

You may want to compare that with a pure niche model (not taking into account the abundances). We can fit the niche models as follow

```
#model comparision
pars_pre_niche <- fit_it(niche_model,
                         Tlevel1 = MPrey, 
                         Tlevel2 = MPred,
                         mean_Tlevel1 = mean(MPrey),
                         sd_Tlevel1 = sd(MPrey),
                         pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                         par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
                         par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10))

plot_pred(pars = pars_pre_niche, Tlevel1 = MPrey, 
          Tlevel2 = MPred, xlab = "log (Predator body size)", ylab = "log (prey body size)",
           pch = ".")
```

We can see in this case that the estimates are similar, but note that the slope (alpha1) is flatter and the range a bit wider. We can compare the accurancy of the integrated, niche and neutral models by comparing its likelyhoods. 

```
#compute likelihoods
lh_model <- -integrated_model(pars_pre, MPrey, MPred, mean(MPrey),sd(MPrey))
lh_niche <- -niche_model(pars_pre_niche, MPrey, MPred, mean(MPrey), sd(MPrey))
lh_neutral <- -neutral_model(MPrey, MPred, mean(MPrey), sd(MPrey))
#visualize it and save it very quick
barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))
l1 <- c("Predation", lh_model, lh_niche, lh_neutral)
```

In this case, the integrated model is a bit better than the niche model, being the neutral model the worst.

Now we do the same for grasshopper hervibore data.

```
# Load data for grasshoppers
grass = read.table("data/Deraison2014.txt", h = TRUE)
head(grass)
#As is experimental data, we subset only the "presences" for ilustrative purposes.
grass <- subset(grass, Herbivory > 0)

#We first fit a binary model, interaction/no interaction.
pars_grass_bin <- fit_it(integrated_model,
                         Tlevel1 = grass$P.Leaf.dry.matter.content, 
                         Tlevel2 = grass$G.Incisive.strength,
                         mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                         sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content)) 
#note that in this case we use the default pars, which are quite unconstrained, but works well in this case.
pars_grass_bin

plot_pred(pars = pars_grass_bin, Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
          Tlevel2 = jitter(grass$G.Incisive.strength), xlab = "Incisive strength", 
          ylab = "Leaf dry matter content")

#compare with the liklyhoods of the three models
pars_grass_bin_niche = fit_it(niche_model, 
                   Tlevel1 = grass$P.Leaf.dry.matter.content, 
                   Tlevel2 = grass$G.Incisive.strength,
                   mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                   sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))
#plot_pred(pars = pars_grass_bin_niche, Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
 #         Tlevel2 = jitter(grass$G.Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")

#likelihoods
lh_model <- -integrated_model(pars_grass_bin, grass$P.Leaf.dry.matter.content, 
                              grass$G.Incisive.strength, 
                              mean(grass$P.Leaf.dry.matter.content), 
                              sd(grass$P.Leaf.dry.matter.content))
lh_niche <- -niche_model(pars_grass_bin_niche, grass$P.Leaf.dry.matter.content, 
                         grass$G.Incisive.strength, 
                         mean(grass$P.Leaf.dry.matter.content), 
                         sd(grass$P.Leaf.dry.matter.content))
lh_neutral <- -neutral_model(grass$P.Leaf.dry.matter.content, grass$G.Incisive.strength, 
                             mean(grass$P.Leaf.dry.matter.content), 
                             sd(grass$P.Leaf.dry.matter.content))

barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))
l2 <- c("Hervibory_bin", lh_model, lh_niche, lh_neutral)
#here the integrated is as good as the neutral model.

#Now we test the frequentistic grasshoper data, taking into account how frequent is each link.
#prepare the data
head(grass)
Incisive.strength <- c()
for(i in 1:nrow(grass)){
  temp <- rep(grass$G.Incisive.strength[i], round(grass$Herbivory[i]))
  Incisive.strength <- c(Incisive.strength,temp)
}
Dry.matter <- c()
for(i in 1:nrow(grass)){
  temp <- rep(grass$P.Leaf.dry.matter.content[i], round(grass$Herbivory[i]))
  Dry.matter <- c(Dry.matter,temp)
}

#fit models
pars_grass_freq <- fit_it(integrated_model,
                         Tlevel1 = Dry.matter, 
                         Tlevel2 = Incisive.strength,
                         mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                         sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))
#note the distribution (mean and sd) is unweigthed because the experiment had equal abundances of plant species. We can use this knowledge to atribute each plant equal weight.
pars_grass_freq
plot_pred(pars = pars_grass_freq, Tlevel1 = jitter(Dry.matter,100), 
          Tlevel2 = jitter(Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")

#compare models
pars_grass_freq_niche <- fit_it(niche_model, 
                                Tlevel1 = Dry.matter, 
                                Tlevel2 = Incisive.strength,
                                mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                                sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))

#plot_pred(pars = pars_grass_freq_niche, Tlevel1 = jitter(Dry.matter,100), 
 #         Tlevel2 = jitter(Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")
#likelihoods
lh_model <- -integrated_model(pars_grass_freq, Dry.matter, Incisive.strength, 
                              mean(grass$P.Leaf.dry.matter.content), 
                              sd(grass$P.Leaf.dry.matter.content))
lh_niche <- -niche_model(pars_grass_freq_niche, Dry.matter, Incisive.strength, 
                         mean(grass$P.Leaf.dry.matter.content), 
                         sd(grass$P.Leaf.dry.matter.content))
lh_neutral <- -neutral_model(Dry.matter, Incisive.strength, 
                             mean(grass$P.Leaf.dry.matter.content),
                             sd(grass$P.Leaf.dry.matter.content))

barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))
l3 <- c("Hervibory_freq", lh_model, lh_niche, lh_neutral)
```

And here is the code for the plants and pollinators:

```
#Read data
pols = read.table("data/Bartomeus2008.txt", h = TRUE)
head(pols)
#Here we transform body size to tongue lenght first based on Cariveau et al. Submitted.
require(devtools)
install_github("BeeIT", "ibartomeus")
library(BeeIT)
head(pols)
pols$tongue <- ITtongue(pols$IT_mm,pols$family , mouthpart = "tongue")
#plot(pols$tongue, pols$IT_mm) #see the relationship if interested
#subset the plant data for later use.
plants <- unique(pols[,c("plant", "nectar_holder_depth_mm", "cover")])

#run models
pars_pols <- fit_it(integrated_model,
                    Tlevel1 = pols$nectar_holder_depth_mm, 
                      Tlevel2 = log(pols$tongue),
                    mean_Tlevel1 = weighted_mean(plants$nectar_holder_depth_mm, plants$cover),
                    sd_Tlevel1 = weighted_sd(plants$nectar_holder_depth_mm, plants$cover))
#note that here we have independent data on plant abundance (% cover). Hence we can use this data directly to fit the model.
#if interested, you can see that the inference from the network, is quite similar to the real cover 
#plot(dnorm(x = seq(-7,10,1), weighted_mean(plants$nectar_holder_depth_mm, plants$cover), 
 #     weighted_sd(plants$nectar_holder_depth_mm, plants$cover)), col = "red", type = "l")
#par(new = TRUE)
#plot(dnorm(x = seq(-7,10,1), mean(pols$nectar_holder_depth_mm), 
 #     sd(pols$nectar_holder_depth_mm)), col = "blue", type = "l")

plot_pred(pars = pars_pols, Tlevel1 = jitter(pols$nectar_holder_depth_mm, 10), 
          Tlevel2 = log(pols$tongue), xlab = "log(Pollinator tongue size)", ylab = "Nectar depth")

#model comparision
pars_pol_niche = fit_it(niche_model, 
                   Tlevel1 = pols$nectar_holder_depth_mm, 
                   Tlevel2 = log(pols$tongue),
                   mean_Tlevel1 = weighted_mean(plants$nectar_holder_depth_mm, plants$cover),
                   sd_Tlevel1 = weighted_sd(plants$nectar_holder_depth_mm, plants$cover))
plot_pred(pars = pars_pol_niche, Tlevel1 = jitter(pols$nectar_holder_depth_mm, 10), 
          Tlevel2 = log(pols$tongue), xlab = "log(Pollinator tongue size)", ylab = "Nectar depth")

#likelihoods
lh_model <- -integrated_model(pars_pols, pols$nectar_holder_depth_mm, log(pols$tongue), 
                              weighted_mean(plants$nectar_holder_depth_mm, plants$cover), 
                              weighted_sd(plants$nectar_holder_depth_mm, plants$cover))
lh_niche <- -niche_model(pars_pol_niche, pols$nectar_holder_depth_mm, log(pols$tongue), 
                         weighted_mean(plants$nectar_holder_depth_mm, plants$cover),
                         sd(pols$nectar_holder_depth_mm))
lh_neutral <- -neutral_model(pols$nectar_holder_depth_mm, log(pols$tongue), 
                             weighted_mean(plants$nectar_holder_depth_mm, plants$cover),
                             weighted_sd(plants$nectar_holder_depth_mm, plants$cover))

barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))
l5 <- c("Pollintion", lh_model, lh_niche, lh_neutral)
```

Hosts and parasitoids. An example that don't work as well.

```
#read and format data
host = read.table("data/Tylianakis2008.txt", h = TRUE)
head(host)
host_body_length <- c()
for(i in 1:nrow(host)){
  temp <- rep(host$host_body_length[i], host$freq[i])
  host_body_length <- c(host_body_length, temp)
}
parasite_body_length <- c()
for(i in 1:nrow(host)){
  temp <- rep(host$parasite_body_length[i], host$freq[i])
  parasite_body_length <- c(parasite_body_length,temp)
}

#models
pars_host <- fit_it(integrated_model,
                        Tlevel1 = host_body_length, 
                        Tlevel2 = parasite_body_length,
                        mean_Tlevel1 = weighted_mean(host$host_body_length, host$freq),
                        sd_Tlevel1 = weighted_sd(host$host_body_length, host$freq),
                        par_lo = c(a0 = 0, a1 = -10, b0 = -10, b1 = -10),
                        par_hi = c(a0 = 10, a1 = 0, b0 = 10, b1 = 10))
#with unconstrained pars, the model behaves wierd, so we have to constrain the values to sensible limits. You can do that by looking at x and y axes ranges and think which slope values are possible.
pars_host
plot_pred(pars = pars_host, Tlevel1 = jitter(host_body_length), 
          Tlevel2 = jitter(parasite_body_length), xlab = "Parasite body size", ylab = "Host body size")
#As you can see everything is predicted.

#compare
pars_host_niche <- fit_it(niche_model,
                    Tlevel1 = host_body_length, 
                    Tlevel2 = parasite_body_length,
                    mean_Tlevel1 = weighted_mean(host$host_body_length, host$freq),
                    sd_Tlevel1 = weighted_sd(host$host_body_length, host$freq),
                    par_lo = c(a0 = 0, a1 = -10, b0 = -10, b1 = -10),
                    par_hi = c(a0 = 10, a1 = 0, b0 = 10, b1 = 10))

plot_pred(pars = pars_host_niche, Tlevel1 = host_body_length, 
          Tlevel2 = parasite_body_length, xlab = "Parasite body size", ylab = "Host body size")

#likelihoods
lh_model <- -integrated_model(pars_host, host_body_length, parasite_body_length, 
                              weighted_mean(host$host_body_length, host$freq), 
                              weighted_sd(host$host_body_length, host$freq))
lh_niche <- -niche_model(pars_host_niche, host_body_length, parasite_body_length, 
                         weighted_mean(host$host_body_length, host$freq), 
                         weighted_sd(host$host_body_length, host$freq))
lh_neutral <- -neutral_model(host_body_length, parasite_body_length, 
                             weighted_mean(host$host_body_length, host$freq), 
                             weighted_sd(host$host_body_length, host$freq))

barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))
l4 <- c("Parasitsim", lh_model, lh_niche, lh_neutral)

```
Thats's it. You can also gather a summary table of parameters and likelyhoods for all models done:

```
#table of likelyhoods
d <- rbind(l1,l2,l3,l4,l5)
d <- as.data.frame(d)
colnames(d) <- c("system", "integrated", "niche", "neutral")
d

dd <- rbind(pars_pre, pars_pre_niche, pars_grass_bin, pars_grass_bin_niche, pars_grass_freq, 
            pars_grass_freq_niche, pars_host, pars_host_niche, pars_pols, pars_pol_niche)
dd <- as.data.frame(dd)
colnames(dd) <- c("a0", "a1", "b0", "b1")
dd

```
