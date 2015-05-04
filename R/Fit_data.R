#Run real data

source('R/FitModel_SA.R')

# Load data for Pred/Prey
#I DIDN'T UPLOAD THE DATA TO GIT BECAUSE THE REPO WILL GO OPEN EVENTUALLY AND I DON'T KNOW IF THIS DATA IS PUBLIC
original_data = read.csv2("temp/size_Barnes2008.csv",dec=".")
head(original_data)
MPred = log10(original_data$standardised_predator_length)
MPrey = log10(original_data$si_prey_length)


#pred/prey----
pars_pre <- fit_it(integrated_model,
                   Tlevel1 = MPred, 
                   Tlevel2 = MPrey,
                   mean_Tlevel1 = mean(MPred),
                   sd_Tlevel1 = sd(MPrey))

plot.pred(pars = pars_pre, Tlevel1 = MPred, 
          Tlevel2 = MPrey, xlab = "log (prey body size)", ylab = "log (Predator body size)")

#model comparision
pars_pre_niche <- fit_it(niche_model,
                   Tlevel1 = MPred, 
                   Tlevel2 = MPrey,
                   mean_Tlevel1 = mean(MPred),
                   sd_Tlevel1 = sd(MPrey))
plot.pred(pars = pars_pre_niche, Tlevel1 = MPred, 
          Tlevel2 = MPrey, xlab = "log (prey body size)", ylab = "log (Predator body size)")

#likelihoods
lh_model <- -integrated_model(pars_pre, MPred, MPrey, mean(MPred),sd(MPrey))
lh_niche <- -niche_model(pars_pre_niche, MPred, MPrey, mean(MPred), sd(MPrey))
lh_neutral <- -neutral_model(MPred, MPrey, mean(MPred), sd(MPrey))

barplot(c(lh_model, lh_niche, lh_neutral))
l1 <- c("Predation", lh_model, lh_niche, lh_neutral)

# Load data for grasshoppers
grass = read.table("data/Deraison2014.txt", h = TRUE)
head(grass)
grass <- subset(grass, Herbivory > 10)

#grasshoppers--------
pars_grass_bin <- fit_it(integrated_model,
                         Tlevel1 = grass$G.Incisive.strength, 
                         Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio,
                         mean_Tlevel1 = mean(grass$G.Incisive.strength),
                         sd_Tlevel1 = sd(grass$G.Incisive.strength))

plot.pred(pars = pars_grass_bin, Tlevel1 = grass$G.Incisive.strength, 
          Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio, xlab = "C:N ratio", ylab = "Incisive strength")
#model comparision
pars_grass_bin_niche = fit_it(niche_model, 
                   Tlevel1 = grass$G.Incisive.strength, 
                   Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio,
                   mean_Tlevel1 = mean(grass$G.Incisive.strength),
                   sd_Tlevel1 = sd(grass$G.Incisive.strength))
plot.pred(pars = pars_grass_bin_niche, Tlevel1 = grass$G.Incisive.strength, 
          Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio, xlab = "C:N ratio", ylab = "Incisive strength")
#likelihoods
lh_model <- -integrated_model(pars_grass_bin, grass$G.Incisive.strength, grass$P.Leaf.carbon.nitrogen.ratio, 
                              mean(grass$G.Incisive.strength), sd(grass$G.Incisive.strength))
lh_niche <- -niche_model(pars_grass_bin_niche, grass$G.Incisive.strength, grass$P.Leaf.carbon.nitrogen.ratio, 
                         mean(grass$G.Incisive.strength), sd(grass$G.Incisive.strength))
lh_neutral <- -neutral_model(grass$G.Incisive.strength, grass$P.Leaf.carbon.nitrogen.ratio, 
                             mean(grass$G.Incisive.strength), sd(grass$G.Incisive.strength))

barplot(c(lh_model, lh_niche, lh_neutral))
l2 <- c("Hervibory_bin", lh_model, lh_niche, lh_neutral)

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
pars_grass_freq <- fit_it(integrated_model, 
                          Tlevel1 = Incisive.strength, 
                          Tlevel2 = Carbon.nitrogen,
                          mean_Tlevel1 = mean(Incisive.strength),
                          sd_Tlevel1 = sd(Incisive.strength))
                          
plot.pred(pars = pars_grass_freq, Tlevel1 = Incisive.strength, 
          Tlevel2 = Carbon.nitrogen, xlab = "C:N ratio", ylab = "Incisive strength")
#Not sure what is wrong!
#check other notation:
pars_grass_freq2 <- fit_it(integrated_model,
                         Tlevel1 = grass$G.Incisive.strength, 
                         Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio,
                         mean_Tlevel1 = weighted.mean(grass$G.Incisive.strength, grass$Herbivory),
                         sd_Tlevel1 = weighted.sd(grass$G.Incisive.strength, grass$Herbivory))

plot.pred(pars = pars_grass_freq2, Tlevel1 = grass$G.Incisive.strength, 
          Tlevel2 = grass$P.Leaf.carbon.nitrogen.ratio, xlab = "C:N ratio", ylab = "Incisive strength")
#bad also, and different :(

pars_grass_freq_niche <- fit_it(niche_model, 
                          Tlevel1 = Incisive.strength, 
                          Tlevel2 = Carbon.nitrogen,
                          mean_Tlevel1 = mean(Incisive.strength),
                          sd_Tlevel1 = sd(Incisive.strength))

plot.pred(pars = pars_grass_freq_niche, Tlevel1 = Incisive.strength, 
          Tlevel2 = Carbon.nitrogen, xlab = "C:N ratio", ylab = "Incisive strength")
#likelihoods
lh_model <- -integrated_model(pars_grass_freq, Incisive.strength, Carbon.nitrogen, 
                              mean(Incisive.strength), sd(Incisive.strength))
lh_niche <- -niche_model(pars_grass_freq_niche, Incisive.strength, Carbon.nitrogen, 
                         mean(Incisive.strength), sd(Incisive.strength))
lh_neutral <- -neutral_model(Incisive.strength, Carbon.nitrogen, 
                             mean(Incisive.strength), sd(Incisive.strength))

barplot(c(lh_model, lh_niche, lh_neutral))
l3 <- c("Hervibory_freq", lh_model, lh_niche, lh_neutral)

#Host/para-------
host = read.table("data/Tylianakis2008.txt", h = TRUE)
head(host)

plot(host$parasite_body_length ~ host$host_body_length)

pars_host <- fit_it(integrated_model,
                        Tlevel1 = host$parasite_body_length, 
                        Tlevel2 = host$host_body_length,
                        mean_Tlevel1 = weighted.mean(host$parasite_body_length, host$freq),
                        sd_Tlevel1 = weighted.sd(host$parasite_body_length, host$freq))

plot.pred(pars = pars_host, Tlevel1 = host$parasite_body_length, 
          Tlevel2 = host$host_body_length, xlab = "Host body size", ylab = "Parasite body size")

#compare models? worth it? the results indicate no trait matching at all.
pars_host_niche <- fit_it(niche_model,
                    Tlevel1 = host$parasite_body_length, 
                    Tlevel2 = host$host_body_length,
                    mean_Tlevel1 = weighted.mean(host$parasite_body_length, host$freq),
                    sd_Tlevel1 = weighted.sd(host$parasite_body_length, host$freq))

plot.pred(pars = pars_host_niche, Tlevel1 = host$parasite_body_length, 
          Tlevel2 = host$host_body_length, xlab = "Host body size", ylab = "Parasite body size")

#likelihoods
lh_model <- -integrated_model(pars_host, host$parasite_body_length, host$host_body_length, 
                              weighted.mean(host$parasite_body_length, host$freq), weighted.sd(host$parasite_body_length, host$freq))
lh_niche <- -niche_model(pars_host_niche, host$parasite_body_length, host$host_body_length, 
                         weighted.mean(host$parasite_body_length, host$freq), weighted.sd(host$parasite_body_length, host$freq))
lh_neutral <- -neutral_model(host$parasite_body_length, host$host_body_length, 
                             weighted.mean(host$parasite_body_length, host$freq), weighted.sd(host$parasite_body_length, host$freq))

barplot(c(lh_model, lh_niche, lh_neutral))
l4 <- c("Parasitsim", lh_model, lh_niche, lh_neutral)


#pollintor----
pols = read.table("data/Bartomeus2008.txt", h = TRUE)
head(pols)
pols[which(pols$nectar_holder_depth_mm > 15),]
pols <- subset(pols, nectar_holder_depth_mm < 15)
#Dominique, note that this outyier has a tremendous influence in the model.
#The outlyer is wrong because Lasioglossums can enter the tube!

#transform to tongue lenght first!

pars_pols <- fit_it(integrated_model,
                    Tlevel1 = pols$IT_mm, 
                    Tlevel2 = pols$nectar_holder_depth_mm,
                    mean_Tlevel1 = mean(pols$IT_mm),
                    sd_Tlevel1 = sd(pols$IT_mm))

plot.pred(pars = pars_pols, Tlevel1 = pols$IT_mm, 
          Tlevel2 = pols$nectar_holder_depth_mm, xlab = "Nectar depth", ylab = "Pollinator size")

#model comparision
pars_pol_niche = fit_it(niche_model, 
                   Tlevel1 = pols$IT_mm, 
                   Tlevel2 = pols$nectar_holder_depth_mm,
                   mean_Tlevel1 = mean(pols$IT_mm),
                   sd_Tlevel1 = sd(pols$IT_mm))
plot.pred(pars = pars_pol_niche, Tlevel1 = pols$IT_mm, 
          Tlevel2 = pols$nectar_holder_depth_mm, xlab = "Nectar depth", ylab = "Pollinator size")
#likelihoods
lh_model <- -integrated_model(pars_pols, pols$IT_mm, pols$nectar_holder_depth_mm, 
                              mean(pols$IT_mm), sd(pols$IT_mm))
lh_niche <- -niche_model(pars_pol_niche, pols$IT_mm, pols$nectar_holder_depth_mm, 
                         mean(pols$IT_mm), sd(pols$IT_mm))
lh_neutral <- -neutral_model(pols$IT_mm, pols$nectar_holder_depth_mm, 
                             mean(pols$IT_mm), sd(pols$IT_mm))

barplot(c(lh_model, lh_niche, lh_neutral))
l5 <- c("Pollintion", lh_model, lh_niche, lh_neutral)

#table of likelyhoods
d <- rbind(l1,l2,l3,l4,l5)
d <- as.data.frame(d)
colnames(d) <- c("system", "integrated", "niche", "neutral")
d
