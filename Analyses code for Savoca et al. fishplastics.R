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
library(lme4)
library(MuMIn)
library(ggrepel)
library(sjPlot)
library(tibble)
library(tidyr)
library(ggeffects)



SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}


# Abbreviate a binomial e.g. Gadus morhua -> G. morhua
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

d_poll <- as_tibble(read_csv("Spatial Information_microplastics.csv"))


d_R1 = read_csv("Plastics ingestion records fish master_final_GCB_v2.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n,
         ProvCode = "oceanographic_province_from_longhurst_2007") %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed"),
         ProvCode = ifelse(ProvCode == "BPRL", "BPLR",
                           ifelse(ProvCode == "HUMB", "CHIL",
                                  ifelse(ProvCode == "NAST E", "NASE",
                                         ifelse(ProvCode == "NPSE", "NPPF", ProvCode))))) %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE) %>% 
  left_join(dplyr::select(d_poll, c(ProvCode, adjacency, mean_poll_abund)), by = "ProvCode") %>% 
  mutate(method_type = factor(method_type),
         adjacency = as_factor(case_when(adjacency == 1 ~ "coastal",
                                         adjacency == 0 ~ "oceanic")),
         source = as_factor(source),
         family = ifelse(family == "Gasterostediae", "Gasterosteidae", 
                         ifelse(family == "Merluccidae", "Merlucciidae", family)),
         adjacency_water = as_factor(ifelse(water_type == "estuarine", "estuarine", 
                                            ifelse(adjacency == "coastal", "coastal", "oceanic")))) %>% 
  
  # Additional edits here for R1
  filter(source != "Sun et al. 2019") %>%   #duplicate with Zhao et al. 2019, removing Sun et al.
  
  mutate(
    equipment_general = case_when(
      equipment_used %in% c("bottom trawls, trammel fishing", "bottom trawls", "bottom trawls, aquaculture", "bottom trawls, purse seine") ~ "bottom trawls",
      equipment_used %in% c("fish market", "fish market, bycatch", "fish market, field collection", "fish market, fishermen") ~ "fish market",
      equipment_used %in% c("fyke nets, gillnets", "gillnets", "gillnets, purse seine", "gillnet, fish market") ~ "gillnets",
      equipment_used %in% c("longline", "line", "nets, hook-and-line", "trolling lines") ~ "lines incld. longlines",
      equipment_used %in% c("trawling", "beam trawl", "fyke nets, trawling", "otter trawl, trawling, longline", "semi-pelagic trawls", "otter trawl", "nets, trawling", "otter trawl, fyke nets", "pelagic trawl", "trawling, bycatch", "pelagic trawl, bottom trawls", "rapido trawl", "shrimp nets") ~ "trawling",
      equipment_used %in% c("purse seine", "beach seine") ~ "seines"),
    
    capture_general = case_when(
      capture_purpose %in% c("Aquaculture", "plastic study, aquaculture") ~ "aquaculture", 
      capture_purpose %in% c("artisanal fishing", "recreational fishing", "artisanal fishing, fish market", "plastic study, artisanal fishery") ~ "small scale fishing",
      capture_purpose %in% c( "bycatch", "incidental") ~ "bycatch",
      capture_purpose %in% c("plastic study (not listed)", "plastic study (not listed)", "plastic study", "plastic study, commercial fishing", "plastic study, fish market", "plastic study, commercial fishing, recreational fishing") ~ "plastic study"),
    
    brk_pt_MD = ifelse(publication_year >= 2017, "after", "before")
  )


d_R1 %>%   filter(prop_w_plastic > 0) %>%  summarise(n_studies = SE(mean_num_particles_per_indv))

#database for all data where microplastics were quantified
d_full <- d_R1 %>%
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  filter(includes_microplastic == "Y") 



# total plastic FO, ALL DATA
sum(d$NwP, na.rm = TRUE)/ sum(d$N, na.rm = TRUE)

SE(d$prop_w_plastic)

# total plastic FO, MP ONLY DATA
sum(d_full$NwP, na.rm = TRUE)/ sum(d_full$N, na.rm = TRUE) 

SE(d_full$prop_w_plastic)

hist(d_full$mean_num_particles_per_indv, xlim = c(0,10), breaks = 10000)

d_full_R1 %>% 
  #filter(average_depth >800) %>% 
  #filter(mean_num_particles_per_indv >0) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = weighted.mean(mean_num_particles_per_indv, w = N, na.rm = TRUE),
            num_plast_SE = SE(mean_num_particles_per_indv),
            med_num_plast = median(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE))



# Expanded dataframe that accounts for sample size with rows----
df_ELH_bp <- d_full_wo_gaps_PF %>% 
  drop_na(N) %>% 
  map_df(., rep, .$N)


PF_EcoMod_bp <- ggplot(df_ELH_bp,
                       aes(x = brk_pt_MD, y = prop_w_plastic)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(height=0.2, width=0.25), alpha = 0.02)
PF_EcoMod_bp



#database for all data where microplastics were quantified; 0 = no, 1 = yes
d_full_R1 <- d_R1 %>%
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  filter(includes_microplastic == "Y") %>% 
  mutate(meth_3_YN = ifelse(method_type != 3, 0,1),
    poly_conf_YN = ifelse(polymer_confirmation %in% c("none", NA), 0,1),
         blanks_used_YN = ifelse(blanks_used %in% c("no", "not described", NA), 0,1),
         clean_lab_YN = ifelse(clean_lab_procedures_used %in% c("no", "not described", NA), 0,1),
         min_size_YN = ifelse(is.na(smallest_detection_size_limits_mm), 0,1),
         overall_reliability = poly_conf_YN + blanks_used_YN + clean_lab_YN + min_size_YN)

d_full_R1_by_study <- d_full_R1 %>% 
  group_by(source) %>% 
  summarise(Overall_FO = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            N = sum(N),
            Method_3 = first(poly_conf_YN),
            Publication_year = first(publication_year),
            Polymer_confirmation = first(poly_conf_YN),
            Blanks_used = first(blanks_used_YN),
            Clean_lab = first(clean_lab_YN),
            Min_size = first(min_size_YN),
            Min_size = replace_na(Min_size, 0),
            Smallest_detect_size = first(smallest_detection_size_limits_mm),
            Method_type = first(method_type)
            # Equipment_general = first(equipment_general),
            # Capture_general = first(capture_general)
  ) %>% 
  mutate(Overall_reliability = Polymer_confirmation + Blanks_used + Clean_lab + Min_size)


gold_standard_studies <- d_full_R1_by_study %>% 
  filter(Polymer_confirmation == 1, Blanks_used == 1, Clean_lab == 1, Min_size == 1,
         Method_type == 3)




# summary tables of raw data----

# species summary table
d_sp_sum <- d_R1 %>%
  filter(!species %in% c("sp.", "spp.","spp")) %>%
  group_by(binomial, family, order, commercial, iucn_status) %>%
  drop_na(binomial, family) %>% 
  summarize(Sp_mean = mean(prop_w_plastic, na.rm = TRUE),
            Sample_size = sum(N, na.rm = TRUE),
            #num_sp = n_distinct(binomial),
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
  # filter(num_sp > 2) %>% 
  filter(Sample_size >10, Sp_mean < 0.25) %>%
  #filter(commercial == "Commercial") %>%
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
conserve_fish <- d_R1 %>% 
  group_by(binomial) %>% 
  #filter(iucn_status == "LC") %>% 
  filter(iucn_status %in% c("NT", "VU", "EN", "CR")) %>% 
  summarize(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            iucn = first(iucn_status),
            Sample_size = sum(N, na.rm = TRUE)) %>% 
  filter(Sample_size > 9) %>% 
  arrange(-FO_plastic)

# Supplementary Table S2 ----
d_vulnerability <- d_full %>%
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
    FO_plastic > 0.25, Sample_size > 10) %>% 
  arrange(-FO_plastic)

write_csv(d_vulnerability, "Vulnerability table.csv")

# Supplementary Table S3, fish of concern for humans ----
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
  filter(species_avg > 0.25 & sample_size >= 10) %>% 
  arrange(-species_avg)
write_csv(concern_fish, "Concerning fish for humans.csv")


# geographic summary of data
Fish_geo_summ <- d_R1 %>% 
  #filter(ProvCode %in% c("CHIN", "KURO", "SUND", "INDE")) %>%
  filter(ProvCode %in% c("NASE", "NECS", "MEDI")) %>% 
  #filter(ProvCode == "BPLR") %>% 
  #filter(water_type != "estuarine") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  #filter(includes_microplastic == "Y") %>% 
  #group_by(ProvCode) %>% 
  summarize(num_studies = n_distinct(source),
            num_sp = n_distinct(binomial),
            num_w_plast = sum(NwP, na.rm = TRUE),
            num_ind_studied = sum(N, na.rm = TRUE),
            prop_by_region = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            wgt_mean_plast_num = weighted.mean(mean_num_particles_per_indv, N),
            se_plast_num = SE(mean_num_particles_per_indv)) 

%>% 
  #left_join(d_poll, by = "ProvCode")

write_csv(Fish_geo_summ, "Fish_plastic geographic summary.csv")

summary(lm(prop_by_region~mean_poll_abund, data = Fish_geo_summ))

# Modeling for fish-plastic ingestion meta-analysis paper----

#GLMM with method type----
FO_MT <- d_full %>% 
  drop_na(N,NwP, order, method_type) 

glmm_FwP_MT <- glmer(cbind(NwP, N-NwP) ~ method_type +
                       (1|order) + (1|source), 
                     na.action = "na.fail",
                     data = FO_MT, family = binomial)
summary(glmm_FwP_MT)
r.squaredGLMM(glmm_FwP_MT) 

# TEST BEHAVIOR, ECOLOGY, GEOGRAPHY----

# full model ----
d_full_wo_gaps_full <- d_full_R1 %>%
  filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  #mutate(prime_forage = ifelse(prime_forage == "Benthic foraging", "aaBenthic foraging", prime_forage)) %>% 
  drop_na(NwP, N, prime_forage, average_depth, Found, trophic_level_via_fishbase) %>% 
  filter(prime_forage != "Scavenging")


glmm_FwP_full <- glmer(cbind(NwP, N-NwP) ~ prime_forage + scale(trophic_level_via_fishbase) + 
                         scale(average_depth)*Found + scale(mean_poll_abund) + adjacency_water +
                         publication_year + brk_pt_MD +
                         (1|order) + (1|method_type) + (1|overall_reliability), 
                     na.action = "na.fail",
                     data = d_full_wo_gaps_full, 
                     family = binomial)
summary(glmm_FwP_full)
r.squaredGLMM(glmm_FwP_full)


GLMM_dredge <- dredge(glmm_FwP_full)

View(glmm_FwP_full)
write_csv(GLMM_dredge, "GLMM model selection table.csv")


d_full_wo_gaps_PF <- d_full_R1 %>%
  filter(includes_microplastic == "Y") %>%
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  mutate(prime_forage = ifelse(prime_forage == "Benthic foraging", "aaBenthic foraging", prime_forage)) %>%
  drop_na(NwP, N, prime_forage) %>%
  filter(prime_forage != "Scavenging")



glmm_FwP_PF_new <- glmer(cbind(NwP, N-NwP) ~ prime_forage + 
                               (1|order)  + (1|publication_year) + (1|method_type) + (1|overall_reliability) + (1|brk_pt_MD),  
                             na.action = "na.fail",
                             data = d_full_wo_gaps_PF, 
                             family = binomial)
summary(glmm_FwP_PF_new)
r.squaredGLMM(glmm_FwP_PF_new)




d_full %>% 
  pull(average_depth) %>% 
  quantile(c(0.025, 0.25, 0.5, 0.75, 0.89, 0.95, 0.99), na.rm = TRUE)

d_full_wo_gaps_AD_F <- d_full_R1 %>%
  filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  #filter(average_depth < 1500) %>%  # Can TOGGLE in and out 
  drop_na(average_depth, Found, NwP, N)

glmm_FwP_AD_F <- glmer(cbind(NwP, N-NwP) ~ 
                         scale(average_depth)*Found + 
                         (1|order)  + (1|publication_year) + (1|method_type) + (1|overall_reliability) + (1|brk_pt_MD), 
                             na.action = "na.fail",
                             data = d_full_wo_gaps_AD_F, 
                             family = binomial)
summary(glmm_FwP_AD_F)
r.squaredGLMM(glmm_FwP_AD_F)

#DO AICc
plot_model(glmm_FwP_AD_F, type = "re",
           grid = FALSE, 
           sort.est = "sort.all")[[2]] +
  ggtitle("") +
  #xlab("Order") +
  ylab("Random intercept") +
  #ylim(0.1,10) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1)

dev.copy2pdf(file="Random effect of order_GLMM.pdf", width=7, height=12)



d_full_wo_gaps_poll <- d_full_R1 %>%
  filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  filter(adjacency_water != "estuarine") %>%
  drop_na(Found, NwP, N, mean_poll_abund)

glmm_FwP_poll  <- glmer(cbind(NwP, N-NwP) ~ 
                          scale(mean_poll_abund) +  #Do separately WITHOUT estuarine
                          (1|order)  + (1|publication_year) + (1|method_type) + (1|overall_reliability) + (1|brk_pt_MD),
                             na.action = "na.fail",
                             data = d_full_wo_gaps_poll, 
                             family = binomial)
summary(glmm_FwP_poll)
r.squaredGLMM(glmm_FwP_poll)



d_full_wo_gaps_AW <- d_full_R1 %>%
  filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  filter(adjacency_water != "estuarine") %>%
  #mutate(adjacency_water = ifelse(adjacency_water == "estuarine", "aaestuarine", adjacency_water)) %>% 
  drop_na(Found, adjacency_water, NwP, N)

glmm_FwP_AW  <- glmer(cbind(NwP, N-NwP) ~ 
                        adjacency_water + 
                        (1|order)  + (1|publication_year) + (1|method_type) + (1|overall_reliability) + (1|brk_pt_MD), 
                        na.action = "na.fail",
                        data = d_full_wo_gaps_AW, 
                        family = binomial)
summary(glmm_FwP_AW)
r.squaredGLMM(glmm_FwP_AW)



#GLMM with trophic level----
d_full_wo_gaps_TL <- d_full_R1 %>%
  #filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  drop_na(NwP, N, trophic_level_via_fishbase)


# This is the ecological model, POSSIBLY LOOK INTO A GAM 
glmm_FwP_TL <- glmer(cbind(NwP, N-NwP) ~ scale(trophic_level_via_fishbase) +
                       (1|order)  + (1|publication_year) + (1|method_type) + (1|overall_reliability) + (1|brk_pt_MD), # AICc with and without pub year; try with FAMILY or ORDER FAMILY AND ORDER CHANGES NUMBERS A LOT
                             na.action = "na.fail",
                             data = d_full_wo_gaps_TL, 
                             family = binomial)
summary(glmm_FwP_TL)
r.squaredGLMM(glmm_FwP_TL)


df_ELH_TL <- d_full_wo_gaps_AD_F %>% 
  drop_na(N) %>% 
  map_df(., rep, .$N)


PF_EcoMod_bp <- ggplot(df_ELH_TL,
                       aes(x = average_depth, y = prop_w_plastic)) +
  geom_smooth() +
  geom_point(alpha = 0.02)
PF_EcoMod_bp




Model_table <- model.sel(glmm_FwP_PF_new, glmm_FwP_AD_F, glmm_FwP_poll, 
                         glmm_FwP_AW, glmm_FwP_TL)
View(Model_table)






# Reviewer 1 test
glmm_FwP_Rev1test <- glmer(cbind(NwP, N-NwP) ~ scale(publication_year) + 
                            (1|brk_pt_MD) + (1|overall_reliability)  + (1|method_type),  
                           na.action = "na.fail",
                           data = d_full_wo_gaps_PF,
                           family = binomial)
summary(glmm_FwP_Rev1test)
r.squaredGLMM(glmm_FwP_Rev2test)





#GLMM with publication year----
FO_year_2010 <- d_full_R1 %>% 
  drop_na(N,NwP, publication_year, order) %>%
  #filter(method_type != 1) %>%  # Can TOGGLE in and out
  #filter(adjacency_water != "estuarine") %>%  # Can TOGGLE in and out
  filter(publication_year >2009)




glmm_FwP_pub_year_R1 <- glmer(cbind(NwP, N-NwP) ~ scale(publication_year) +
                                (1|brk_pt_MD) + (1|overall_reliability)  + (1|order) + (1|method_type), 
                           na.action = "na.fail",
                           data = FO_year_2010, family = binomial)
summary(glmm_FwP_pub_year_R1)
r.squaredGLMM(glmm_FwP_pub_year_R1) 

FO_year_lm <- lm(prop_w_plastic~publication_year, weights = N, data = FO_year_2010)
summary(FO_year_lm)

