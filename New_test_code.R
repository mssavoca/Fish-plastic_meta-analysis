# Trying revisions data (updated w four comments from reviewers 1)

library(tidyverse)
library(ggpubr)



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
         

d_R1 %>% summarise(n_studies = n_distinct(source))



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
    mutate(poly_conf_YN = ifelse(polymer_confirmation %in% c("none", NA), 0,1),
           blanks_used_YN = ifelse(blanks_used %in% c("no", "not described", NA), 0,1),
           clean_lab_YN = ifelse(clean_lab_procedures_used %in% c("no", "not described", NA), 0,1),
           min_size_YN = ifelse(is.na(smallest_detection_size_limits_mm), 0,1),
           overall_reliability = poly_conf_YN + blanks_used_YN + clean_lab_YN + min_size_YN)

d_full_R1_by_study <- d_full_R1 %>% 
  group_by(source) %>% 
  summarise(Overall_FO = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            N = sum(N),
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
  filter(Clean_lab == 1, Blanks_used == 1, Polymer_confirmation ==1, Method_type == 3)



# Exploratory plots for new data ----

Equipment_bp <- ggplot(d_full_R1, aes(equipment_general, prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  labs(x = "Equipement used",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw(base_size = 18)
Equipment_bp


Capture_bp <- ggplot(d_full_R1, aes(capture_general, prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  labs(x = "Capture purpose",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw()
Capture_bp


PolyConf_bp <- ggplot(d_full_R1, aes(as.factor(poly_conf_YN), prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  scale_x_discrete(labels = c('No','Yes')) +
  labs(x = "Polymer confirmation",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw(base_size = 18)
PolyConf_bp


Blanks_bp <- ggplot(d_full_R1, aes(as.factor(blanks_used_YN), prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  scale_x_discrete(labels = c('No','Yes')) +
  labs(x = "Blank methods",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw(base_size = 18)
Blanks_bp



Clean_bp <- ggplot(d_full_R1, aes(as.factor(clean_lab_YN), prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  scale_x_discrete(labels = c('No','Yes')) +
  labs(x = "Clean lab procedures described",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw(base_size = 18)
Clean_bp


MinSize_bp <- ggplot(d_full_R1, aes(as.factor(min_size_YN), prop_w_plastic)) +
  geom_boxplot() +
  geom_point(aes(size = N),
             position=position_jitter(height=0, width=0.2), alpha = 0.5) +
  scale_x_discrete(labels = c('No','Yes')) +
  labs(x = "Minimum detection (size) threshold reported",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_bw(base_size = 18)
MinSize_bp

  
#trying binomial plots for study quality----
  
bp_polyconf <- ggplot(filter(d_full_R1_by_study, Publication_year >2009),
              aes(Publication_year, Polymer_confirmation)) +
  stat_smooth(method="glm", method.args=list(family="binomial"),formula=y~x, 
              alpha=0.2, size=1) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  labs(x = "",
       y = "Pr(polymer confirmation)") +
  scale_x_continuous(breaks = 2010:2020) +
  theme_minimal(base_size = 12)
bp_polyconf


bp_blanks <- ggplot(filter(d_full_R1_by_study, Publication_year >2009),
                      aes(Publication_year, Blanks_used)) +
  stat_smooth(method="glm", method.args=list(family="binomial"),formula=y~x, 
              alpha=0.2, size=1) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  labs(x = "",
       y = "Pr(blanks used)") +
  scale_x_continuous(breaks = 2010:2020) +
  theme_minimal(base_size = 12)
bp_blanks


bp_cleanlab <- ggplot(filter(d_full_R1_by_study, Publication_year >2009),
                    aes(Publication_year, Clean_lab)) +
  stat_smooth(method="glm", method.args=list(family="binomial"),formula=y~x, 
              alpha=0.2, size=1) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  labs(x = "Publication year",
       y = "Pr(clean lab procedures)") +
  scale_x_continuous(breaks = 2010:2020) +
  theme_minimal(base_size = 12)
bp_cleanlab


bp_reportmin <- ggplot(filter(d_full_R1_by_study, Publication_year >2009),
                      aes(Publication_year, Min_size)) +
  stat_smooth(method="glm", method.args=list(family="binomial"),formula=y~x, 
              alpha=0.2, size=1) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  labs(x = "Publication year",
       y = "Pr(report minimum detection threshold)") +
  scale_x_continuous(breaks = 2010:2020) +
  theme_minimal(base_size = 12)
bp_reportmin




# Combined binomial plots----


BP_Fig_full <- ggarrange(bp_polyconf, bp_blanks, 
                         bp_cleanlab, bp_reportmin,
                        labels = c("A","B","C","D"), 
                        font.label = list(size = 14),
                        align = "hv",
                        legend = "none",
                        ncol = 2, nrow = 2)
BP_Fig_full





bp_all <- ggplot(filter(d_full_R1_by_study, Publication_year >2009)) +
  
  
  # stat_smooth(aes(Publication_year, Method_3),
  #             color = "aquamarine3",
  #             method="glm", method.args=list(family="binomial"),formula=y~x,
  #             alpha=0.2, size=1) +
  # geom_point(aes(Publication_year, Method_3),
  #            position=position_jitter(height=0.05, width=0.05),
  #            color = "aquamarine3", alpha = 0.2) +
  
  stat_smooth(aes(Publication_year, Min_size),
              color = "black",
              method="glm", method.args=list(family="binomial"),formula=y~x,
              alpha=0.2, size=1) +
  geom_point(aes(Publication_year, Min_size),
             position=position_jitter(height=0.05, width=0.05), 
             color = "black", alpha = 0.2) +
  
  stat_smooth(aes(Publication_year, Clean_lab),
              method="glm", method.args=list(family="binomial"),formula=y~x, 
              color = "dodgerblue3",alpha=0.2, size=1) +
  geom_point(aes(Publication_year, Clean_lab),
             position=position_jitter(height=0.05, width=0.05),
             color  = "dodgerblue3", alpha = 0.2) +
  
  stat_smooth(aes(Publication_year, Blanks_used),
              method="glm", method.args=list(family="binomial"),formula=y~x,
              color = "firebrick2",alpha=0.2, size=1) +
  geom_point(aes(Publication_year, Blanks_used),
             position=position_jitter(height=0.05, width=0.05),
             color  = "firebrick2", alpha = 0.2) +

  stat_smooth(aes(Publication_year, Polymer_confirmation),
              method="glm", method.args=list(family="binomial"),formula=y~x,
              color = "chocolate2",alpha=0.2, size=1) +
  geom_point(aes(Publication_year, Polymer_confirmation),
             position=position_jitter(height=0.05, width=0.05),
             color  = "chocolate2", alpha = 0.2) +
  
  #geom_vline(xintercept = 2015, color = "gray40", linetype = "dashed") +
  
  labs(x = "Publication year",
       y = "Quality assurance \nmetrics described") +
  scale_x_continuous(breaks = 2010:2020) +
  scale_y_continuous(breaks = c(0,1), labels = c("No","Yes")) +
  theme_classic(base_size = 14) 
bp_all



# color palette for figures
pal <- c("Min_size" = "black", "Clean_lab" = "dodgerblue3", "Blanks_used" = "firebrick2",  "Polymer_confirmation" = "chocolate2")


bp_all_long <- d_full_R1_by_study %>% 
 filter(Publication_year >2009) %>%  #can toggle in and out
  pivot_longer(cols = Polymer_confirmation:Min_size, names_to = "Metric_type", values_to = "Metric_value") %>% 

  ggplot() + 
  
  stat_smooth(aes(Publication_year, Metric_value, color = Metric_type),
              method="glm", method.args=list(family="binomial"),formula=y~x, 
              alpha=0.2, size=1) +
  geom_point(aes(Publication_year, Metric_value, color = Metric_type),
             position=position_jitter(height=0.05, width=0.05),
             alpha = 0.2) +
  
  #geom_vline(xintercept = 2015, color = "gray40", linetype = "dashed") +
  
  scale_color_manual(values = pal, labels = c("Blanks \nused", 
                                              "Clean lab methods \ndescribed", 
                                              "Minimum size \nthreshold reported", 
                                              "Polymer \nconfirmation")) +
  
  labs(x = "Publication year",
       y = "Quality assurance \nmetrics described",
       color = "Quality assurance \nmetric") +
  scale_x_continuous(breaks = 2010:2020) +
  scale_y_continuous(breaks = c(0,1), labels = c("No","Yes")) +
  theme_classic(base_size = 12) +
  theme(legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(0.75, "cm"),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

bp_all_long

dev.copy2pdf(file="Quality_assurance_plot.pdf", width=5.25, height=4)


#plot of minimum size of plastic----

Detect_FO <- d_full_R1_by_study %>% 
  mutate(Publication_year = as_factor(Publication_year)) %>% 
  ggplot() +
  geom_point(
    aes(Smallest_detect_size, Overall_FO,
                     weight = N, size= N, color = Publication_year), 
    alpha = 0.5) +
  scale_colour_viridis_d(option = "viridis") +

  geom_smooth(data = filter(d_full_R1_by_study, Smallest_detect_size > 0.15),
              aes(Smallest_detect_size, Overall_FO, weight = N, size= N),
              method = "lm") +

  geom_smooth(data = filter(d_full_R1_by_study, Smallest_detect_size < 0.25),
              aes(Smallest_detect_size, Overall_FO, weight = N, size= N),
              formula = y ~ x + I(x^2),
              level = 0.90,
              method = "lm") +
  
  # geom_smooth(data = d_full_R1_by_study,
  #             aes(Smallest_detect_size, Overall_FO),
  #   method="gam", formula = y ~ x + I(x^2), se = F) +
  # geom_hline(yintercept = 0.26, linetype="dashed", color = "grey50") +
  # scale_color_gradientn(colours = c("steelblue4",
  #                                   "darkgoldenrod1",
  #                                   "darkorange", "orangered1",
  #                                   "firebrick1", "red3", "red4")) +
  scale_size_continuous(breaks = c(1, 10, 100, 500, 1000)) +
  scale_x_reverse(breaks=seq(0, 1, 0.05), limits = c(1, 0)) +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "black") +
  #scale_shape_manual(na.translate = F, values = shapes) +
  ylim(-0.1, 1) +
  labs(x = "Smallest size detected (mm)",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
Detect_FO

dev.copy2pdf(file="FO_by_min_size.pdf", width=7.5, height=5)



Detect_PubYear <- ggplot(filter(d_full_R1_by_study, Publication_year >2009),
                         aes(Publication_year, Smallest_detect_size, 
                             weight = N, size= N)) +
  scale_x_reverse() +
  geom_point() +
  geom_smooth(method = "lm") +

  scale_size_continuous(breaks = c(1, 10, 100, 500, 1000)) +
  scale_x_continuous(breaks=seq(2010, 2020, 1)) +
  #scale_shape_manual(na.translate = F, values = shapes) +
  ylim(0, 1) +
  labs(x = "Publication year",
       y = "Smallest size detected (mm)",
       size = "Sample size") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
Detect_PubYear






# Two y-axes plot ----

d_full_R1_tap <- d_full_R1 %>% 
  mutate(smallest_detection_size_limits_mm = rev(smallest_detection_size_limits_mm))


Detect_FO_PubYear <- ggplot() +

  geom_point(data = filter(d_full_R1, publication_year >2009),
             aes(publication_year, prop_w_plastic, 
             alpha = 0.2, size = N, weight = N), color = "indianred1") +
  geom_smooth(data = filter(d_full_R1, publication_year >2009),
              aes(publication_year, prop_w_plastic,
                  size = N, weight = N),
              method = "lm", color = "red3") +
  
  geom_point(data = filter(d_full_R1_by_study, Publication_year >2009),
             aes(Publication_year, 1-Smallest_detect_size), 
             alpha = 0.5, shape = 18, color = "blue") +
  geom_smooth(data = filter(d_full_R1_by_study, Publication_year >2009),
              aes(Publication_year, 1-Smallest_detect_size),
              method = "lm", color = "blue") +
  
  #scale_size_continuous(breaks = c(1, 10, 100, 500, 1000)) +
  scale_x_continuous(breaks=seq(2010, 2020, 1)) +
  scale_y_continuous(
    name = "Proportion with ingested plastic",
    sec.axis = sec_axis(~rev(.), name = "Minimum size thresold (mm)")) +
  
  labs(x = "Publication year",
       size = "Sample size") +
  theme_classic(base_size = 16) +
  theme(
   # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y.left = element_text(color = "red3"),
        axis.text.y.left = element_text(color = "red3"),
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none")
Detect_FO_PubYear


dev.copy2pdf(file="Detect_FO_PubYear_tap.pdf", width=7, height=5)

 # Trying 3D plot ----
Try_3D <- plot_ly(
  x = d_full_R1_by_study$Publication_year, 
  y = d_full_R1_by_study$Overall_FO, 
 # z = d_full_R1_by_study$Smallest_detect_size, 
  #size = N,
) %>%
  #add_markers() %>%
  layout(
    scene = list(xaxis = list(title = 'Publication year'),
                 yaxis = list(title = 'Proportion with ingested plastic'),
                 zaxis = list(title = 'Smallest size detected (mm)'))
  ) %>% 
  add_surface(z = volcano)
Try_3D



Try_3D <- rgl::persp3d(x = d_full_R1_by_study$Publication_year,
                  y = d_full_R1_by_study$Overall_FO, 
                  z = d_full_R1_by_study$Smallest_detect_size, 
                  col="skyblue")
Try_3D



# Figure 1, ORDER phylogeny----

d_order <- d_full %>% 
  group_by(order) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_sp = n_distinct(binomial),
            num_studies = n_distinct(source),
            prop_commercial = sum(commercial %in% c("commercial", "highly commercial"))/n()) %>% 
  drop_na(order) %>% 
  mutate(studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial_cat = cut(prop_commercial,
                              breaks=c(-Inf, 0.01, 0.25, Inf), 
                              labels=c("None", "Minor", "Commercial"))) %>% 
  filter(FO_plastic >0.25, sample_size >9, num_sp > 2) %>% 
  arrange(desc(FO_plastic))


d_order <- d_full %>% 
  group_by(order) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_studies = n_distinct(source),
            prop_commercial = sum(commercial %in% c("commercial", "highly commercial"))/n()) %>% 
  drop_na(order) %>% 
  mutate(studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial_cat = cut(prop_commercial,
                              breaks=c(-Inf, 0.01, 0.25, Inf), 
                              labels=c("None", "Minor", "Commercial"))) %>% 
  arrange(desc(FO_plastic))

# adding
taxon_search <- tnrs_match_names(names = d_order$order, context_name = "All life")
ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
node_in_tree <- node_id(taxon_search)[is_in_tree(node_id(taxon_search))]
d_order$ott_name <- unique_name(taxon_search)
d_order$ott_id <- taxon_search$ott_id



#plots ORDERS
my_tree <- tol_induced_subtree(ott_ids = c(d_order$ott_id),
                               label_format = "name")

my_tree$tip.label[my_tree$tip.label == "mrcaott17080ott19944"] <- "Carcharhiniformes"
my_tree$tip.label[my_tree$tip.label == "Pristiformes/Rhiniformes_group"] <- "Squaliformes"
my_tree$tip.label[my_tree$tip.label == "Carangaria"] <- "Perciformes"
my_tree$tip.label[my_tree$tip.label == "Gasterosteales"] <- "Trachiniformes"
my_tree$tip.label[my_tree$tip.label == "mrcaott13841ott64982"] <- "Anguilliformes"



# LOOK INTO adding Scombriformes for Tuna !!!

a = as_tibble(my_tree)

my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

#add.species.to.genus(my_tree, "Alepisauridae_Synodontidae", where = "root")

# adding familes to taxa manually 

node.tree<-function(tree, m = 0, prefix = "NodE"){
  res <- list(call = match.call())
  res$m.start <- m
  if(is.null(tree$node.label)){
    tree <- ape::makeNodeLabel(tree)
    tree$node.label <- rep("", tree$Nnode)
  }
  for(i in 1:tree$Nnode){
    if(tree$node.label[i] == "NA" | tree$node.label[i] == ""){
      m <- m+1
      tree$node.label[i] <- paste(prefix, m, sep = "")
    }
  }
  res$m.current <- m
  res$prefix <- prefix
  res$tree <- tree
  return(res)
}
tree.label.info <- function(tree, label){
  if(is.null(tree$edge.length) | any(is.na(tree$edge.length))){
    stop("\n tree$edge.length with any NA or NULL\n")
  }
  if(is.null(tree$tip.label) | is.null(tree$node.label)){
    stop("\n tree$tip.label and/or tree$node.label NULL\n")
  }
  n <- length(tree$tip.label)
  where <- which(c(tree$tip.label, tree$node.label) == label)
  if(length(where)>1){
    stop("\n Only one label accepted or label with multiple occurrences in tree\n")
  }
  if(length(where)!=1){
    stop("\n label not found\n")
  }
  H <- phytools::nodeHeights(tree)
  HM <- max(H)
  res <- data.frame(row.names = label)
  if(where<=n){
    res$edge <- which(tree$edge[,2] == where)
    res$edge.length <- tree$edge.length[res$edge]
    res$edge.height <- HM - H[res$edge, 2]
    res$max.height <- HM
    res$type <- "tip"
  }else{
    if (where == (length(tree$tip.label) + 1)){
      res$edge<-which(tree$edge[,1] == where)[1]
      res$edge.length <- 0
      res$edge.height <- HM - H[res$edge, 1]
      res$edge <- NA
      res$max.height <- HM
      res$type <- "root"
    } else {
      res$edge <- which(tree$edge[,2] == where)
      res$edge.length <-tree$edge.length[res$edge]
      res$edge.height <- HM - H[res$edge, 2]
      res$max.height <- HM
      res$type <- "node"
    }
  }
  return(res)
}
add.taxa.phylo <- function(tree, taxa, m = 0, prefix = "NeWNodEPhylO"){
  res <- list(call = match.call())
  res$m.start <- m
  if (!inherits(tree, "phylo")){ 
    stop("\n Object tree is not of class phylo \n")
  }
  if(is.null(tree$node.label)){
    stop("\n tree$node.label is NULL. Use the function node.tree\n")
  }
  tree.labels <- c(tree$tip.label, tree$node.label)
  match.names <- match(taxa[,1], tree.labels)
  if(any(is.na(match.names))){
    print("Some nodes/tips in taxa are not present in tree$node.label/tree$tip.label:")
    mNA <- is.na(match.names)
    for(i in (1:nrow(taxa))[mNA]){
      print(as.character(taxa[i, 1]))
    }
    stop("\n Check tree and/or taxa")
  }
  show.warning <- FALSE
  for(i in 1:nrow(taxa)){
    if(length(which(c(tree.labels) == taxa[i, 1]))>1){
      stop(paste("\n The label", taxa[i,1], "with multiple occurrences in tree$node.label"))
    }
    info.temp <- tree.label.info(tree, taxa[i, 1])
    if(!is.na(taxa[i,3]) & (info.temp$type == "node" | info.temp$type == "root")){
      taxa[i,3] <- NA
      show.warning <- TRUE
    }
    if(!is.na(taxa[i,3]) & info.temp$edge.length<taxa[i,3]){
      stop(paste("\n The edge length to larger to tip", taxa[i,2]))
    }
  }
  if(show.warning){
    warning("\n Edge length not used in internal node anchor")
  }
  u.taxa <- unique(taxa[,1])
  for(i in 1:length(u.taxa)){
    if(length(unique(as.character(taxa[which(taxa[,1] == u.taxa[i]),3])))>1)
      stop("\n Edge length in a sigle anchor tip must be equal")
  }
  edges.length <- as.numeric(taxa[, 3, drop = FALSE])
  control<-matrix(NA, nrow(taxa), 2)
  control[,1] <- taxa[,1]
  for(i in 1:nrow(taxa)){
    label <- taxa[i, 1]
    if(length(which(label == tree$tip.label)) > 0){
      control.temp<-control[which(control[,1] == label), 2]
      if(all(is.na(control.temp))){
        where <- which(tree$tip.label == label)
        if(is.na(edges.length[i])){
          edges.length[i] <- tree.label.info(tree, label)[1, 2]/2 #edge.length
        }
        tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = edges.length[i])
      } else {
        control.tips <- c(control.temp, label)
        control.tips <- control.tips[!is.na(control.tips)]
        node.name <- c(tree$tip.label, tree$node.label)[ape::getMRCA(tree, control.tips)]
        where <- which(c(tree$tip.label, tree$node.label) == node.name)
        tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = 0)
      }
      control[i, 2] <- taxa[i, 2]
    } else {
      where <- which(c(tree$tip.label, tree$node.label) == label)
      tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = 0)
    }
    tree.temp <- node.tree(tree, m = m, prefix = prefix)
    tree <- tree.temp$tree
    m <- tree.temp$m.current
  }
  res$m.current <- m
  res$prefix <- prefix
  res$new.tips <- taxa[,2]
  res$tree <- tree
  return(res)
}


# 
# my_tree <- add.taxa.phylo(my_tree, taxa)$tree
# my_tree <- as.phylo(my_tree) 
# 
# a = as_tibble(my_tree)
# 
#THINGS TO ADD IN: Pleuronectiformes, Scorpaeniformes, Gasterosteiformes (Should be where "Trachiniformes" now is?)

final_taxa <- matrix(c("Zeiformes", "Tetraodontiformes", "Perciformes",  # This row exists on phylogeny
                       "Gasterosteiformes", "Pleuronectiformes", "Scorpaeniformes",   # These are tips to add on phylogeny
                       NA, NA, NA),
                     3,3)
final_taxa

my_tree <- add.taxa.phylo(my_tree, final_taxa)$tree
my_tree <- as.phylo(my_tree)



my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

# first plot try, fan layout
p <- ggtree(my_tree, layout="fan", open.angle=160) +
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  geom_tiplab2(parse = TRUE,
               size = 9,offset = 0.05) +
  ggplot2::xlim(-0.6, 1.3) 


p <- rotate_tree(p, -12)
p

#dev.copy2pdf(file="test_phylo.pdf", width=20, height=20)

# Adding data
shapes <- c("None" = 15, "Minor" = 17, "Commercial" = 16)

p %<+% d_order + 
  aes(color = FO_plastic) +
  geom_tippoint(aes(color = FO_plastic, shape = commercial_cat, size = studies_cat)) +
  scale_color_gradientn(colours = c("steelblue4", "darkgoldenrod1", 
                                    "darkorange", "orangered2", "red3", "red4"), 
                        name = "Proportion with \ningested plastic") +
  #scale_size(range = c(3, 7)) +
  #scale_size_continuous(guide = FALSE, range = c(3, 7)) +
  scale_size_continuous(breaks = seq(from = 1, to = 3, by = 1), 
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n>3)"),
                        range = c(3, 10)) +
  scale_shape_manual(na.translate = F, values = shapes) +
  labs(shape = "Commercial \nstatus") +
  theme(legend.position = c(0.49, 0.49),
        legend.key.size = unit(1.25, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.box = "horizontal") +
  guides(
    size = FALSE, shape = FALSE, color  = FALSE
    #shape = guide_legend(override.aes = list(size = 5))
  )



dev.copy2pdf(file="Fish_order_plastic_phylo_final_nolegend_NEE.pdf", width=22, height=20)
ggsave("Fish_order_plastic_phylo_final_nolegend_NEE.jpg", width = 20, height = 20, units = "in")










MarkicStudies <- read_csv("Markic2019_studymaster.csv") %>% 
  dplyr::select(source, `Method type`) %>% 
  arrange(source)


%>% 
  left_join(d, by = "source")



d_full_wo_gaps_AW_MPPI <- d %>%
  #filter(includes_microplastic == "Y") %>% 
  filter(method_type == 3) %>%  # Can TOGGLE in and out
  #filter(adjacency_water != "estuarine") %>%
  mutate(adjacency_water = ifelse(adjacency_water == "estuarine", "aaestuarine", adjacency_water)) %>% 
  drop_na(mean_num_particles_per_indv, order, Found, adjacency_water, NwP, N)

glmm_FwP_AW_MPPI  <- lmer(log(mean_num_particles_per_indv) ~ 
                            adjacency_water +  #Maybe do separately WITHOUT estuarine
                            (1|order) + (1|source),  # METHOD TYPE ADDED TO INCLUDE AS MUCH DATA AS POSSIBLE
                          na.action = "na.fail",
                          data = d_full_wo_gaps_AW_MPPI)
summary(glmm_FwP_AW_MPPI)
r.squaredGLMM(glmm_FwP_AW_MPPI)



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



# wrapper fucntion  needed to allow MuMIn to compare models
gamm4 <- function(...) structure(c(gamm4::gamm4(...), list(call = match.call())), class = c("gamm", "list"))  

# GAM(M) models
gamm4_FwP_w_conditional <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + 
                                   s((average_depth), by = adjacency_water) + 
                                   Found, 
                                 random = ~(1|order) + (1|source), 
                                 data = filter(d_full_wo_gaps_TL, average_depth < 1501), 
                                 family = binomial)

summary(gamm4_FwP_w_conditional$gam)
plot(gamm4_FwP_w_conditional$gam)



gam_FwP_w_conditional <- gam(prop_w_plastic ~ s(trophic_level_via_fishbase) +
                               s(average_depth, by = Found) +
                               adjacency_water,
                             weights = N,
                             data = filter(d_full_wo_gaps_TL, average_depth < 301),  
                             method = "REML")
summary(gam_FwP_w_conditional)
plot(gam_FwP_w_conditional, all.terms = TRUE, pages = 1)

# GAM model sel
model.sel(gamm4_FwP_w_conditional, gam_FwP_w_conditional, glmm_FwP_eco_geo_TL)


model.sel(gamm4_FwP_w_conditional, gam_FwP_w_conditional, glmm_FwP_eco_geo_TL)

# multi-model selection using AICc  **THINK ABOUT MAKING DEPTH PLOT**
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
  filter(includes_microplastic == "Y", method_type != 1) %>% 
  drop_na(NwP, N, adjacency_water, mean_poll_abund, order, source) 

glmm_FwP_geo <- glmer(cbind(NwP, N-NwP) ~ scale(mean_poll_abund) + 
                        adjacency_water +
                        (1|order) + (1|source), 
                      na.action = "na.fail",
                      data = d_full_wo_gaps_geo, family = binomial)
summary(glmm_FwP_geo)
r.squaredGLMM(glmm_FwP_geo) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  

model.sel(glmm_FwP_eco_geo_PF,glmm_FwP_eco_geo_TL, glmm_FwP_geo)


#for BRTs


###BRT code----
library(dismo)
library(gbm)

try(BRT_FwP_PF <- gbm.step(data=as.data.frame(d_full_BRT),
                           gbm.x=c("prime_forage","order", "average_depth", "Found"), 
                           gbm.y="prop2", 
                           #distribution = "bernoulli",
                           #gbm.y=cbind("NwP", "N"-"NwP"), 
                           #family="binomial",
                           tree.complexity=5, 
                           learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))
gbm.plot(BRT_FwP_PF, write.title = FALSE)
#quartz(); gbm.plot.fits(BRT_FwP_PF)
summary(BRT_FwP_PF)
