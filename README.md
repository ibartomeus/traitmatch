# trait_match
Package to predict trait matching from species interactions.


#Run real data

source('R/FitModel.R')

# Load data for Pred/Prey
#I DIDN'T UPLOAD THE DATA TO GIT BECAUSE THE REPO WILL GO OPEN EVENTUALLY AND I DON'T KNOW IF THIS DATA IS PUBLIC
original_data = read.csv2("temp/size_Barnes2008.csv",dec=".")
head(original_data)
#MPred = log10(original_data$standardised_predator_length)
#MPrey = log10(original_data$si_prey_length)
#cause problems with the precision of the integral of the normal distribution. 
#It has very strange behaviour when the standard deviation is small (<<1). 
#solved by *10

MPred = log10(original_data$standardised_predator_length*10)
MPrey = log10(original_data$si_prey_length*10)

#pred/prey----
pars_pre <- fit_it(integrated_model,
                   Tlevel1 = MPrey, 
                   Tlevel2 = MPred,
                   mean_Tlevel1 = mean(MPrey),
                   sd_Tlevel1 = sd(MPrey))

plot.pred(pars = pars_pre, Tlevel1 = MPrey, 
          Tlevel2 = MPred, xlab = "log (Predator body size)", ylab = "log (prey body size)")

#model comparision
pars_pre_niche <- fit_it(niche_model,
                         Tlevel1 = MPrey, 
                         Tlevel2 = MPred,
                         mean_Tlevel1 = mean(MPrey),
                         sd_Tlevel1 = sd(MPrey))

plot.pred(pars = pars_pre_niche, Tlevel1 = MPrey, 
          Tlevel2 = MPred, xlab = "log (Predator body size)", ylab = "log (prey body size)")

#likelihoods
lh_model <- -integrated_model(pars_pre, MPrey, MPred, mean(MPrey),sd(MPrey))
lh_niche <- -niche_model(pars_pre_niche, MPrey, MPred, mean(MPrey), sd(MPrey))
lh_neutral <- -neutral_model(MPrey, MPred, mean(MPrey), sd(MPrey))

barplot(c(lh_model, lh_niche, lh_neutral))
l1 <- c("Predation", lh_model, lh_niche, lh_neutral)

# Load data for grasshoppers
grass = read.table("data/Deraison2014.txt", h = TRUE)
head(grass)
grass <- subset(grass, Herbivory > 0)

plot(grass$P.Leaf.carbon.nitrogen.ratio, grass$G.Incisive.strength)
plot(grass$P.Leaf.dry.matter.content ~ jitter(grass$G.Incisive.strength), pc = 19,  col = as.numeric(cut(grass$Herbivory, breaks = 4)))

#grasshoppers--------
pars_grass_bin <- fit_it(integrated_model,
                         Tlevel1 = grass$P.Leaf.dry.matter.content, 
                         Tlevel2 = grass$G.Incisive.strength,
                         mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                         sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))

plot.pred(pars = pars_grass_bin, Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
          Tlevel2 = jitter(grass$G.Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")
#model comparision
pars_grass_bin_niche = fit_it(niche_model, 
                   Tlevel1 = grass$P.Leaf.dry.matter.content, 
                   Tlevel2 = grass$G.Incisive.strength,
                   mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                   sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))
plot.pred(pars = pars_grass_bin_niche, Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
          Tlevel2 = jitter(grass$G.Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")
#likelihoods
lh_model <- -integrated_model(pars_grass_bin, grass$P.Leaf.dry.matter.content, grass$G.Incisive.strength, 
                              mean(grass$P.Leaf.dry.matter.content), sd(grass$P.Leaf.dry.matter.content))
lh_niche <- -niche_model(pars_grass_bin_niche, grass$P.Leaf.dry.matter.content, grass$G.Incisive.strength, 
                         mean(grass$P.Leaf.dry.matter.content), sd(grass$P.Leaf.dry.matter.content))
lh_neutral <- -neutral_model(grass$P.Leaf.dry.matter.content, grass$G.Incisive.strength, 
                             mean(grass$P.Leaf.dry.matter.content), sd(grass$P.Leaf.dry.matter.content))

barplot(c(lh_model, lh_niche, lh_neutral))
l2 <- c("Hervibory_bin", lh_model, lh_niche, lh_neutral)

# Make frequentistic grasshoper data
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

plot.pred(pars = pars_grass_freq, Tlevel1 = jitter(Dry.matter,100), 
          Tlevel2 = jitter(Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")
#compare
pars_grass_freq_niche <- fit_it(niche_model, 
                                Tlevel1 = Dry.matter, 
                                Tlevel2 = Incisive.strength,
                                mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                                sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content))

plot.pred(pars = pars_grass_freq_niche, Tlevel1 = jitter(Dry.matter,100), 
          Tlevel2 = jitter(Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")
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

barplot(c(lh_model, lh_niche, lh_neutral))
l3 <- c("Hervibory_freq", lh_model, lh_niche, lh_neutral)

#Host/para-------
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

#plot(host$parasite_body_length ~ host$host_body_length)

pars_host <- fit_it(integrated_model,
                        Tlevel1 = host_body_length, 
                        Tlevel2 = parasite_body_length,
                        mean_Tlevel1 = weighted.mean(host$host_body_length, host$freq),
                        sd_Tlevel1 = weighted.sd(host$host_body_length, host$freq))
#check
weighted.mean(host$host_body_length, host$freq) == mean(host_body_length)

plot.pred(pars = pars_host, Tlevel1 = jitter(host_body_length), 
          Tlevel2 = jitter(parasite_body_length), xlab = "Parasite body size", ylab = "Host body size")

#compare
pars_host_niche <- fit_it(niche_model,
                    Tlevel1 = host_body_length, 
                    Tlevel2 = parasite_body_length,
                    mean_Tlevel1 = weighted.mean(host$host_body_length, host$freq),
                    sd_Tlevel1 = weighted.sd(host$host_body_length, host$freq))

plot.pred(pars = pars_host_niche, Tlevel1 = host_body_length, 
          Tlevel2 = parasite_body_length, xlab = "Parasite body size", ylab = "Host body size")

#likelihoods
lh_model <- -integrated_model(pars_host, host_body_length, parasite_body_length, 
                              weighted.mean(host$host_body_length, host$freq), weighted.sd(host$host_body_length, host$freq))
lh_niche <- -niche_model(pars_host_niche, host_body_length, parasite_body_length, 
                         weighted.mean(host$host_body_length, host$freq), weighted.sd(host$host_body_length, host$freq))
lh_neutral <- -neutral_model(host_body_length, parasite_body_length, 
                             weighted.mean(host$host_body_length, host$freq), weighted.sd(host$host_body_length, host$freq))

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
require(devtools)
install_github("BeeIT", "ibartomeus")
library(BeeIT)
head(pols)
pols$family <- as.character(pols$pollinator)
pols$family[which(pols$family == "Halictus_sp")] <- "Halictidae"
pols$family[which(pols$family == "Nomada_sp")] <- "Apidae"
pols$family[which(pols$family == "Lasioglossum_sp")] <- "Halictidae"
pols$family[which(pols$family == "Andrena_sp")] <- "Andrenidae"
pols$family[which(pols$family == "Panurgus_sp")] <- "Andrenidae"
pols$family[which(pols$family == "Anthophora_sp")] <- "Apidae"
pols$family[which(pols$family == "Osmia_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Apis_mellifera")] <- "Apidae"
pols$family[which(pols$family == "Colletes_sp")] <- "Colletidae"
pols$family[which(pols$family == "Hoplitis??_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Halictus_pyrenaicus")] <- "Halictidae"
pols$family[which(pols$family == "Andrena_agilissima")] <- "Andrenidae"
pols$family[which(pols$family == "Panurgus_banksianus")] <- "Andrenidae"
pols$family[which(pols$family == "Andrena_cinerea")] <- "Andrenidae"
pols$family[which(pols$family == "Prosopis(Hylaeus)_sp??")] <- "Colletidae"
pols$family[which(pols$family == "Lasioglossum??_sp")] <- "Halictidae"
pols$family[which(pols$family == "Lasioglossum_smaragdulus??")] <- "Halictidae"
pols$family[which(pols$family == "Bombus_terrestris")] <- "Apidae"
pols$family[which(pols$family == "Andrena_ovatula")] <- "Andrenidae"
pols$family[which(pols$family == "Hoplitis_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Ceratina_cucurbitina")] <- "Megachilidae"
pols$family[which(pols$family == "Anthophora_plumipes")] <- "Apidae"
pols$family[which(pols$family == "Anthidium_septemspinosum")] <- "Megachilidae"
pols$family[which(pols$family == "Lasioglossum_gemmeus")] <- "Halictidae"
pols$family[which(pols$family == "Megachile_pilicrus")] <- "Megachilidae"
pols$family[which(pols$family == "Anthidium_sticticum")] <- "Megachilidae"
pols$family[which(pols$family == "Osmia_cornuta")] <- "Megachilidae"
pols$family[which(pols$family == "Halictus??_sp")] <- "Halictidae"
pols$family[which(pols$family == "Andrena_nigrocianea??")] <- "Andrenidae"
pols$family[which(pols$family == "Osmia_aurulenta")] <- "Megachilidae"
pols$family[which(pols$family == "Anthidium_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Nomada_succincta")] <- "Apidae"
pols$family[which(pols$family == "Hylaeus_nigritus")] <- "Halictidae"
pols$family[which(pols$family == "Eucera_sp")] <- "Apidae"
pols$family[which(pols$family == "Andrena_carbonaria?")] <- "Andrenidae"
pols$family[which(pols$family == "Xylocopa_sp")] <- "Apidae"
pols$family[which(pols$family == "Bombus_sp")] <- "Apidae"
pols$family[which(pols$family == "Chalicodoma_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Hylaeus_variegatus")] <- "Halictidae"
pols$family[which(pols$family == "Hylaeus_bisinuatus")] <- "Halictidae"
pols$family[which(pols$family == "Eucera_collaris")] <- "Apidae"
pols$family[which(pols$family == "Eumenes_sp")] <- "Apidae"
pols$family[which(pols$family == "Prosopis(Hylaeus)_variegatus")] <- "Halictidae"
pols$family[which(pols$family == "Andrena??_sp")] <- "Andrenidae"
pols$family[which(pols$family == "Dasypoda_cingulata")] <- "Andrenidae"
pols$family[which(pols$family == "Xylocopa_iris")] <- "Apidae"
pols$family[which(pols$family == "Epeolus_variegatus")] <- "Apidae"
pols$family[which(pols$family == "Hylaeus_nigritum")] <- "Halictidae"
pols$family[which(pols$family == "Megachile_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Halictus_gemmeus")] <- "Halictidae"
pols$family[which(pols$family == "-")] <- "Apidae"
pols$family[which(pols$family == "Philanthus_pulchellus")] <- "Apidae"
pols$family[which(pols$family == "Ceratina_dentiventris")] <- "Megachilidae"
pols$family[which(pols$family == "Megachile??_sp")] <- "Megachilidae"
pols$family[which(pols$family == "Halictus_quadricinctus")] <- "Halictidae"
pols$family[which(pols$family == "Prosopis(Hylaeus)_sp")] <- "Colletidae"
pols$family[which(pols$family == "Halictus_pyrenaicus??")] <- "Halictidae"
pols$family[which(pols$family == "Anthidium_stictium")] <- "Megachilidae"
pols$family[which(pols$family == "Halictus_smaragdulus??")] <- "Halictidae"
pols$family[which(pols$family == "Sphecodes_sp")] <- "Halicatidae"
pols$family[which(pols$family == "Hylaeus_sp")] <- "Colletidae"
pols$family[which(pols$family == "Halicatidae")] <- "Halictidae"
pols$family[which(pols$family == "_")] <- "Apidae"

which(!pols$family %in% c('Andrenidae', 'Apidae', 
                         'Colletidae', 'Halictidae', 'Megachilidae'))

pols$tongue <- ITtongue(pols$IT_mm,pols$family , mouthpart = "tongue")

plot(pols$tongue, pols$IT_mm)

#run models
pars_pols <- fit_it(integrated_model,
                    Tlevel1 = pols$nectar_holder_depth_mm, 
                      Tlevel2 = log(pols$tongue),
                    mean_Tlevel1 = mean(pols$nectar_holder_depth_mm),
                    sd_Tlevel1 = sd(pols$nectar_holder_depth_mm))

plot.pred(pars = pars_pols, Tlevel1 = jitter(pols$nectar_holder_depth_mm, 10), 
          Tlevel2 = log(pols$tongue), xlab = "log(Pollinator tongue size)", ylab = "Nectar depth")

#model comparision
pars_pol_niche = fit_it(niche_model, 
                   Tlevel1 = pols$nectar_holder_depth_mm, 
                   Tlevel2 = log(pols$tongue),
                   mean_Tlevel1 = mean(pols$nectar_holder_depth_mm),
                   sd_Tlevel1 = sd(pols$nectar_holder_depth_mm))
plot.pred(pars = pars_pol_niche, Tlevel1 = jitter(pols$nectar_holder_depth_mm, 10), 
          Tlevel2 = log(pols$tongue), xlab = "log(Pollinator tongue size)", ylab = "Nectar depth")
#likelihoods
lh_model <- -integrated_model(pars_pols, pols$nectar_holder_depth_mm, log(pols$tongue), 
                              mean(pols$nectar_holder_depth_mm), sd(pols$nectar_holder_depth_mm))
lh_niche <- -niche_model(pars_pol_niche, pols$nectar_holder_depth_mm, log(pols$tongue), 
                         mean(pols$nectar_holder_depth_mm), sd(pols$nectar_holder_depth_mm))
lh_neutral <- -neutral_model(pols$nectar_holder_depth_mm, log(pols$tongue), 
                             mean(pols$nectar_holder_depth_mm), sd(pols$nectar_holder_depth_mm))

barplot(c(lh_model, lh_niche, lh_neutral))
l5 <- c("Pollintion", lh_model, lh_niche, lh_neutral)

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

