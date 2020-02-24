# Script to replicate analysis for Savoca et al., fish-plastic ingestion meta-analysis paper----

#load packages, data, and functions ----

library(tidyverse)
library(ggtree)
library(ggstance)
library(ape)
library(phylobase)
library(dismo)
library(readxl)
library(rotl) 
library(phytools)
library(tidyverse)
library(stringr)
library(mgcv)
library(gamm4)
library(MuMIn)
library(ggrepel)
library(sjPlot)
library(tibble)
library(tidyr)
library(ggeffects)

SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}


# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

d_poll <- as_tibble(read_csv("Spatial Information_microplastics.csv"))

d = read_csv("Plastics ingestion records fish master_final.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n) %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed")) %>% 
  rename(ProvCode = "oceanographic_province_from_longhurst_2007") %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE) %>% 
  left_join(dplyr::select(d_poll, c(ProvCode, adjacency, mean_poll_abund)), by = "ProvCode") %>% 
  mutate(adjacency = as_factor(case_when(adjacency == 1 ~ "coastal",
                                         adjacency == 0 ~ "oceanic")),
         source = as_factor(source),
         family = ifelse(family == "Gasterostediae", "Gasterosteidae", 
                         ifelse(family == "Merluccidae", "Merlucciidae", family))) 

# total plastic FO
sum(d$NwP, na.rm = TRUE)/ sum(d$N, na.rm = TRUE)            

#database for all data where microplastics were quantified
d_full <- d %>%
  filter(includes_microplastic == "Y") 

# total plastic FO
sum(d_full$NwP, na.rm = TRUE)/ sum(d_full$N, na.rm = TRUE) 
# summary tables of raw data----

# species summary table
d_sp_sum <- d %>%
  filter(!species %in% c("sp.", "spp.","spp")) %>%
  group_by(binomial, family, order, commercial, iucn_status) %>%
  drop_na(binomial, family) %>% 
  summarize(Sp_mean = mean(prop_w_plastic, na.rm = TRUE),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>% 
  ungroup %>% 
  mutate(commercial = factor(commercial),
         studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial = fct_collapse(commercial,
                                   Commercial = c("commercial", "highly commercial"),
                                   Minor = c("minor commercial", "subsistence"),
                                   None = "none")) %>%
  arrange(-Sp_mean)



d_family_disc <- d_full %>% 
  filter(family %in%  c("Carangidae", "Mugilidae", "Pleuronectidae", "Soleidae", "Myctophidae")) %>% 
  group_by(family) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_sp = n_distinct(binomial), 
            num_studies = n_distinct(source)) %>% 
  arrange(desc(FO_plastic))
write_csv(d_family_disc, "Fish families of concern.csv")

# IUCN and vulnerability status of fish
conserve_fish <- d %>% 
  group_by(binomial) %>% 
  #filter(iucn_status == "LC") %>% 
  filter(iucn_status %in% c("NT", "VU", "EN", "CR")) %>% 
  summarize(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            iucn = first(iucn_status),
            Sample_size = sum(N, na.rm = TRUE)) %>% 
  filter(Sample_size > 25) %>% 
  arrange(-FO_plastic)

# Supplementary Table S2
d_vulnerability <- d_full %>%
  #  filter(family %in% c("Scombridae", "Sphyrnidae", "Carcharhinidae")) %>% 
  drop_na(vulnerability_score_via_fishbase_from_cheug_et_al_2005) %>% 
  mutate(Vulnerability.category = cut(vulnerability_score_via_fishbase_from_cheug_et_al_2005,
                                      breaks=c(-Inf, 20, 40, 60, 80, Inf), 
                                      labels=c("low","low-moderate", "moderate-high", "high-very high", "very high"))) %>% 
  group_by(binomial, common_name, vulnerability_score_via_fishbase_from_cheug_et_al_2005, iucn_status) %>%  
  summarize(FO_plastic = sum(NwP)/sum(N),
            Sample_size = sum(N),
            #    num_sp = n_distinct(binomial), 
            num_studies = n_distinct(source)) %>% 
  ungroup %>% 
  filter(
    vulnerability_score_via_fishbase_from_cheug_et_al_2005 >50, 
    #Vulnerability.category %in% c("moderate-high", "high-very high", "very high") & 
    #iucn_status %in% c("NT","VU", "EN", "CR") & 
    FO_plastic > 0.25, Sample_size > 25, ) %>% 
  arrange(-FO_plastic)

write_csv(d_vulnerability, "Vulnerability table.csv")

# Supplementary Table S3, fish of concern for humans
concern_fish <- d %>% 
  group_by(common_name, binomial, family) %>% 
  filter(commercial %in% c("commercial", "highly commercial")) %>%
  summarize(species_avg = sum(NwP)/sum(N),
            mean_num_plast = weighted.mean(mean_num_particles_per_indv, N),
            sample_size = sum(N),
            num_studies = n_distinct(source),
            commercial_status = first(commercial),
            AC_status = first(aquaculture), 
            rec_status = first(recreational)) %>% 
  ungroup %>% 
  filter(species_avg > 0.25 & sample_size > 25) %>% 
  arrange(-species_avg)
write_csv(concern_fish, "Concerning fish for humans.csv")


# geographic summary of data
Fish_geo_summ <- d_full %>% 
  filter(ProvCode %in% c("CHIN", "KURO", "SUND", "INDE")) %>%
  #filter(ProvCode %in% c("NAST E", "NECS", "MEDI")) %>% 
  #filter(ProvCode == "BPRL") %>% 
  group_by(ProvCode) %>% 
  summarize(num_studies = n_distinct(source),
            num_sp = n_distinct(binomial),
            num_w_plast = sum(NwP, na.rm = TRUE),
            num_ind_studied = sum(N, na.rm = TRUE),
            prop_by_region = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            wgt_mean_plast_num = weighted.mean(mean_num_particles_per_indv, N),
            se_plast_num = SE(mean_num_particles_per_indv))
write_csv(Fish_geo_summ, "Fish_plastic geographic summary.csv")



# Modeling for fish-plastic ingestion meta-analysis paper----

#GLMM with publication year----
FO_year_2010 <- d_full %>% 
  drop_na(N,prop_w_plastic, publication_year, order) %>% 
  filter(publication_year >2009) 

glmm_FwP_pub_year <- glmer(cbind(NwP, N-NwP) ~ scale(publication_year) +
                             (1|order) + (1|source), 
                           na.action = "na.fail",
                           data = FO_year_2010, family = binomial)
summary(glmm_FwP_pub_year)
r.squaredGLMM(glmm_FwP_pub_year) 

# Testing which full ecological/geographic model is better---
model.sel(glmm_FwP_eco_geo_TL, glmm_FwP_eco_geo_PF)

#Full GLMM with foraging behavior, which is better model by AIC----

d_full_wo_gaps_PF <- d %>%
  filter(includes_microplastic == "Y") %>% 
  mutate(prime_forage = ifelse(prime_forage == "Benthic foraging", "aaBenthic foraging", prime_forage)) %>% 
  drop_na(average_depth, Found, NwP, N, prime_forage,
          adjacency, mean_poll_abund, ProvCode)

glmm_FwP_eco_geo_PF <- glmer(cbind(NwP, N-NwP) ~ scale(average_depth) + 
                               scale(mean_poll_abund) +
                               prime_forage + Found + adjacency +
                               (1|order) + (1|source), 
                             na.action = "na.fail",
                             data = d_full_wo_gaps_PF, family = binomial)
summary(glmm_FwP_eco_geo_PF)
r.squaredGLMM(glmm_FwP_eco_geo_PF)

# multi-model selection using AICc
GLMM_dredge_PF <- dredge(glmm_FwP_eco_geo_PF)

View(GLMM_dredge_PF)
write_csv(GLMM_dredge_PF, "GLMM model selection table.csv")
PF_avg = model.avg(GLMM_dredge_PF)
summary(PF_avg) #The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. An alternative, the ‘full’ average assumes that a variable is included in every model


# plot random effects WOW
plot_model(glmm_FwP_eco_geo_PF, type = "re",
           grid = FALSE, 
           sort.est = "sort.all")[[1]] +
  ggtitle("") +
  xlab("Reference") +
  ylab("Random intercept") +
  #ylim(0.1,10) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  geom_hline(yintercept = 1)

dev.copy2pdf(file="Random effect of source_GLMM.pdf", width=7, height=11)

plot_model(glmm_FwP_eco_geo_PF, type = "re",
           grid = FALSE, 
           sort.est = "sort.all")[[2]] +
  ggtitle("") +
  xlab("Order") +
  ylab("Random intercept") +
  #ylim(0.1,10) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1)

dev.copy2pdf(file="Random effect of order_GLMM.pdf", width=7, height=11)


# MCMCglmm_FwP_eco_geo_PF <- MCMCglmm(prop_w_plastic ~  prime_forage + average_depth + mean_poll_abund + Found + adjacency,
#                                     random = ~source,
#                                     data = d_full_wo_gaps_PF,
#                                     family = "gaussian",
#                                     nitt = 11000, thin = 1, burnin = 1000,
#                                     pr = TRUE, # To save the posterior mode of for each level or category in your random effect(s) we use pr = TRUE, which saves them in the $Sol part of the model output.
#                                     verbose = TRUE)
# summary(MCMCglmm_FwP_eco_geo_PF)
# 
# MCMCglmm_FwP_eco_geo_PF$Fixed


#Full GLMM with trophic level----
d_full_wo_gaps_TL <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, trophic_level_via_fishbase, NwP, N, 
          adjacency, mean_poll_abund, ProvCode) 



glmm_FwP_eco_geo_TL <- glmer(cbind(NwP, N-NwP) ~ scale(trophic_level_via_fishbase) +
                               scale(average_depth) +
                               scale(mean_poll_abund) +
                               Found + 
                               adjacency + 
                               (1|order) + (1|source), 
                             na.action = "na.fail",
                             data = d_full_wo_gaps_TL, family = binomial)
summary(glmm_FwP_eco_geo_TL)
r.squaredGLMM(glmm_FwP_eco_geo_TL)




# multi-model selection using AICc
GLMM_dredge <- dredge(glmm_FwP_eco_geo_TL)

View(GLMM_dredge)
write_csv(GLMM_dredge, "GLMM model selection table.csv")
#subset(GLMM_dredge, delta < 4)
a=model.avg(GLMM_dredge)
summary(a) #The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. An alternative, the ‘full’ average assumes that a variable is included in every model
confint(a) #computes confidence interval


# Comparing Eco-only to Geo-only models----
glmm_FwP_eco <- glmer(cbind(NwP, N-NwP) ~ prime_forage + 
                        scale(average_depth) + 
                        Found + 
                        (1|order) + (1|source), 
                      na.action = "na.fail",
                      data = d_full_wo_gaps_PF, family = binomial)
summary(glmm_FwP_eco)
r.squaredGLMM(glmm_FwP_eco) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  


d_full_wo_gaps_geo <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(NwP, N, adjacency, mean_poll_abund, order, source) 

glmm_FwP_geo <- glmer(cbind(NwP, N-NwP) ~ scale(mean_poll_abund) + adjacency +
                        (1|order) + (1|source), 
                      na.action = "na.fail",
                      data = d_full_wo_gaps_geo, family = binomial)
summary(glmm_FwP_geo)
r.squaredGLMM(glmm_FwP_geo) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  

model.sel(glmm_FwP_eco, glmm_FwP_geo)
