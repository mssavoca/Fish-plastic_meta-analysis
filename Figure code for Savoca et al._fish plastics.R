# Script to replicate figures for Savoca et al., fish-plastic ingestion paper----


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
library(ggpubr)

SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}


# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

d_poll <- as_tibble(read_csv("Spatial Information_microplastics.csv"))


d = read_csv("Plastics ingestion records fish master_final_SciAd.csv") %>% 
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
                                            ifelse(adjacency == "coastal", "coastal", "oceanic"))))


#database for all data where microplastics were quantified
d_full <- d %>%
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  filter(includes_microplastic == "Y") 




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


#Code for figures ----

# Figure 1, Family phylogeny----

d_family <- d_full %>% 
  group_by(family) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_sp = n_distinct(binomial),
            num_studies = n_distinct(source),
            prop_commercial = sum(commercial %in% c("commercial", "highly commercial"))/n()) %>% 
  drop_na(family) %>% 
  mutate(studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial_cat = cut(prop_commercial,
                              breaks=c(-Inf, 0.01, 0.25, Inf), 
                              labels=c("None", "Minor", "Commercial"))) %>% 
  #filter(FO_plastic >0.25, sample_size >25, num_sp > 2, commercial_cat == "Commercial") %>% 
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
taxon_search <- tnrs_match_names(names = d_family$family, context_name = "All life")
ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
node_in_tree <- node_id(taxon_search)[is_in_tree(node_id(taxon_search))]
d_family$ott_name <- unique_name(taxon_search)
d_family$ott_id <- taxon_search$ott_id

#plots FAMILIES
my_tree <- tol_induced_subtree(ott_ids = c(85219, 118776, 816156, 19305, 893503, 960231, 563241, 253371, 486695, 42415, 738324, 561376, 239745, 679814, 749780, 34184, 875698, 
                                           978560, 1003121, 384676, 137651, 739933, 191025, 646019, 804461, 722754, 1089734, 637234, 648754, 37461, 888446, 563518, 856584,
                                           186486, 563230, 479853, 734459, 437596, 1089730, 441564, 292707, 823203, 427507, 42408, 655424, 44805, 577720, 19301, 441571, 232621, 
                                           279762, 1089742, 813991, 1032209, 129124, 804451, 308750, 587923, 583638, 340592, 400235, 749770, 1074732, 701559, 563513, 710014, 
                                           99942, 450143, 609233, 160291, 712841, 739939, 883406, 1089740, 460871, 734202, 818997, 214115, 132684, 72393, 407171, 563254, 74014, 
                                           280947, 99937, 769569, 214115, 175436, 308741, 114163, 363181, 615333, 258647, 1052881, 13838, 555245, 1093612, 930712, 778687, 579429, 
                                           562630, 574822, 280953, 339090, 540474, 859881, 137656, 62639, 811925, 479864, 644001, 95055, 401063, 765787, 715629, 724438, 659851, 
                                           614519, 769567, 892951, 892958, 237354, 250743, 65336, 553102, 1026498, 36225),
                               label_format = "name")


my_tree$tip.label[my_tree$tip.label == "Belonidae"] <- "Scomberesocidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott801ott480916"] <- "Ariidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott595ott3629079"] <- "Sternoptychidae"
my_tree$tip.label[my_tree$tip.label == "Gonostomatidae_(family_in_Opisthokonta)"] <- "Gonostomatidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott47507ott47509"] <- "Acanthuridae" 
my_tree$tip.label[my_tree$tip.label == "mrcaott7529ott77068"] <- "Cottidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott9658ott106653"] <- "Pholidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott15200ott303035"] <- "Nototheniidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott126948ott301378"] <- "Cheilodactylidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott7815ott9211"] <- "Gobiidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott61966ott80330"] <- "Bythitidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott56885ott4134216"] <- "Lotidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott13841ott182974"] <- "Nettastomatidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott65747ott167810"] <- "Ophichthidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott73302ott93287"] <- "Atherinopsidae"
my_tree$tip.label[my_tree$tip.label == "Alepocephaliformes"] <- "Alepocephalidae"
my_tree$tip.label[my_tree$tip.label == "Ephippiformes"] <- "Ephippidae"

a = as_tibble(my_tree)
which(my_tree$tip.label=="Lutjanidae")

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

taxa <- matrix(c("Squalidae", "Squalidae", "Sphyrnidae", "mrcaott17080ott24097", "mrcaott17080ott24097", "Alepisauridae", "Lotidae", "Engraulidae",  "Muraenidae", "Sebastidae", "Pleuronectidae", "mrcaott7429ott39883", "Alepisauridae", "Gobiidae", "Sillaginidae", "Terapontidae", "mrcaott595ott7526", "Moronidae", "Priacanthidae", "Nototheniidae", "Cottales", "Trichiuridae",  "Moridae", "mrcaott497ott21417", "Trichiuridae", "Alepisauridae", "Exocoetidae", "mrcaott3549ott7508", "Lateolabracidae", "mrcaott595ott7526",  # This row exists on phylogeny
                 "Somniosidae", "Etmopteridae", "Carcharhinidae", "Scyliorhinidae", "Triakidae", "Synodontidae", "Gadidae", "Clupeidae","Muraenesocidae","Scorpaenidae", "Paralichthyidae", "Acropomatidae", "Chlorophthalmidae", "Eleotridae", "Sparidae", "Kyphosidae", "Stomiidae", "Serranidae", "Cepolidae", "Eleginopsidae", "Hexagrammidae", "Scombridae", "Macrouridae", "Carangidae", "Gempylidae", "Paralepidae", "Hemiramphidae", "Scaridae", "Polyprionidae", "Phosichthyidae",   # These are tips to add on phylogeny
                 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), 
               30,3)
taxa

my_tree <- add.taxa.phylo(my_tree, taxa)$tree
my_tree <- as.phylo(my_tree) 

a = as_tibble(my_tree)

#THINGS TO ADD IN IN SECOND PHYLO CALL
#Congridae # sister to Muraenesocidae # ADD THIS IN AT THE END IN A SEPARATE ADD.TAXA CALL
# Labridae # sister to Scaridae # ADD THIS IN AT THE END IN A SEPARATE ADD.TAXA CALL

final_taxa <- matrix(c("Muraenesocidae", "Scaridae", "mrcaott33200ott33222",  # This row exists on phylogeny
                       "Congridae", "Labridae", "Merlucciidae",   # These are tips to add on phylogeny
                       NA,NA,NA), 
                     3,3)
final_taxa

my_tree <- add.taxa.phylo(my_tree, final_taxa)$tree
my_tree <- as.phylo(my_tree) 



my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

# first plot try, fan layout
p <- ggtree(my_tree, layout="circular", open.angle=90) +
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  geom_tiplab2(parse = TRUE,
               size = 6,offset = 0.05) +
  ggplot2::xlim(-0.6, 1.3) 
p

#dev.copy2pdf(file="test_phylo.pdf", width=20, height=20)

# Adding data
shapes <- c("None" = 15, "Minor" = 17, "Commercial" = 16)

p %<+% d_family + 
  aes(color = FO_plastic) +
  # geom_tiplab2(aes(label = paste0("italic('", label, "')"),
  #                  color=FO_plastic), parse = TRUE,
  #              size = 6, align = FALSE, offset = 0.05) +
  geom_tippoint(aes(color = FO_plastic, shape = commercial_cat, size = studies_cat)) +
  scale_color_gradientn(colours = c("steelblue4", "darkgoldenrod1", 
                                    "darkorange", "orangered2", "red3", "red4"), 
                        name = "Proportion with \ningested plastic") +
  #scale_size(range = c(3, 7)) +
  #scale_size_continuous(guide = FALSE, range = c(3, 7)) +
  scale_size_continuous(breaks = seq(from = 1, to = 3, by = 1), 
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n>3)"),
                        range = c(3, 7)) +
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

dev.copy2pdf(file="Fish_family_plastic_phylo_final_nolegend_SciAd.pdf", width=20, height=20)
ggsave("Fish_family_plastic_phylo_final_nolegend_SciAd.jpg", width = 20, height = 20, units = "in")


# Figure 2 (code also for S2a,b and S3) ----
##Map depicting frequency of occurrence of plastic ingestion according to Longhurst province. 

#Here we also show data preparation for most maps in this publication (Figures 2, S3 - S5). This code bins and organizes the data for use in QGIS. The resulting .csv was exported from R and imported into QGIS.  More detailed instructions on how to import .csv files into QGIS are listed at the end of this section. 

# load packages and data
library(tidyverse)
library(gbm)
library(dismo)
library(mgcv)
library(lme4)
library(gamm4)
library(readxl)
library(readr)
library(ggplot2)

### For all figures (2, S2a and S2b, S3)
dat <- read.csv("/Users/test/Box Sync/Microplastics/Fish-plastic_meta-analysis/Plastics ingestion records fish master_final_SciAd.csv")
# get rid of estuarine studies
data <- dat[dat$Water.type=="marine",] # Use this full data set for Figures S2a and S3
nrow(dat)
nrow(data)

## Figure S2b: include this line if you're looking for studies that use Method 3
data <- data %>% filter(Method.type==3) 
###################

### Figure 2 ###
# include this line if you're looking for studies that have to include microplastics
data <- data %>% filter(Includes.microplastic == "Y"|Includes.microplastic == "Y?") 



# Here we want to get average proportion of plastic per province
# data of interest is:
data2 <- data[,c("Binomial", "Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "NwP", "N", "Source")]
colnames(data2) <- c("Species", "OceanProv", "PropPlastic", "NwP", "N", "Source")
data3<- data2[order(data2$OceanProv),]

# We need our provinces to match the publicly available shape file, so we are renaming those that don't and then dropping the replaced factor levels
levels(data3$OceanProv) <- c(levels(data3$OceanProv), "CHIL", "BPLR", "NASE", "NPPF") 
data3$OceanProv[data3$OceanProv=="BPRL"] <- "BPLR"
data3$OceanProv[data3$OceanProv=="HUMB"] <- "CHIL"
data3$OceanProv[data3$OceanProv=="NAST E"] <- "NASE"
data3$OceanProv[data3$OceanProv=="NPSE"] <- "NPPF"

# some checks to make sure that worked
table(data3$OceanProv)
data3$OceanProv <- droplevels(data3$OceanProv)
levels(data3$OceanProv)

# Here we will use a for-loop to fill in the desired variables (average plastic, number of fish, number of studies, number of species, normalized sample size per province)

#Set up new dataframe to include our desired variables
prov <- unique(data3$OceanProv)
prov
length(prov)
aveplast <- rep(NA, length(prov))
numfish <- rep(NA, length(prov))
numstudies <- rep(NA, length(prov))
numspecies <- rep(NA, length(prov))
normalized <- rep(NA, length(prov))

# Now our for-loop
for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- sum(sub$NwP, na.rm=T)/sum(sub$N, na.rm=T)
  numfish[i] <- sum(sub$N, na.rm=T)
  numstudies[i] <- length(unique(sub$Source))
  numspecies[i] <- length(unique(sub$Species))
  normalized[i] <- aveplast[i]*sum(sub$N, na.rm=T)
}

# Checking output
aveplast
length(aveplast)
numfish
table(numfish)
numstudies
head(data3)
numspecies

#Combine all of these variables into a new data frame, which we will build off of to make our map
newdat <- data.frame(prov, aveplast, numfish, numstudies, numspecies, normalized)
newdat

# We are going to bin data for easier categorization and comparison:
# ave plastic
newdat$aveplastbin <- rep(NA, length(prov))
newdat$aveplastbin[newdat$aveplast<=.10]= "<.10"
newdat$aveplastbin[newdat$aveplast >.10 & newdat$aveplast<=20]= ".11-.20"
newdat$aveplastbin[newdat$aveplast >.20 & newdat$aveplast<=.30]= ".21-.30"
newdat$aveplastbin[newdat$aveplast >.30 & newdat$aveplast<=.40]= ".31-.40"
newdat$aveplastbin[newdat$aveplast >.40 & newdat$aveplast<=.50]= ".41-.50"
newdat$aveplastbin[newdat$aveplast >.50 & newdat$aveplast<=.60]= ".51-.60"
newdat$aveplastbin[newdat$aveplast >.60 & newdat$aveplast<=.70]= ".61-.70"
newdat$aveplastbin[newdat$aveplast >.70 & newdat$aveplast<=.80]= ".71-.80"
newdat$aveplastbin[newdat$aveplast >.80 & newdat$aveplast<=.90]= ".81-.90"
newdat$aveplastbin[newdat$aveplast >.90 & newdat$aveplast<=.99]= ".91-.99"
newdat$aveplastbin[newdat$aveplast ==1]= "1"

# number of studies
summary(newdat$numstudies)
newdat$numstudiesbin <- rep(NA, length(prov))
newdat$numstudiesbin[newdat$numstudies == 1] = "1"
newdat$numstudiesbin[newdat$numstudies >1 & newdat$numstudies<=5]= "2-5"
newdat$numstudiesbin[newdat$numstudies >5 & newdat$numstudies<=10]= "6-10"
newdat$numstudiesbin[newdat$numstudies >=10]= ">10"

# number of fish/study
summary(newdat$numfish)
newdat$numfishbin <- rep(NA, length(prov))
newdat$numfishbin[newdat$numfish >=1 & newdat$numfish<10]= "< 10"
newdat$numfishbin[newdat$numfish >=10 & newdat$numfish<=50]= "10-50"
newdat$numfishbin[newdat$numfish >50 & newdat$numfish<=100]= "51-100"
newdat$numfishbin[newdat$numfish >100 & newdat$numfish<500]= "101-500"
newdat$numfishbin[newdat$numfish >500 & newdat$numfish<=1000]= "501-1000"
newdat$numfishbin[newdat$numfish >1000 & newdat$numfish<=1500]= "1001-1500"
newdat$numfishbin[newdat$numfish >1500]= ">1500"

# normalized proportion of ingestion bins
summary(newdat$normalized)
newdat$normbin <- rep(NA, length(prov))
newdat$normbin[newdat$normalized >=1 & newdat$normalized<25]= "< 10"
newdat$normbin[newdat$normalized >=10 & newdat$normalized<=50]= "10-50"
newdat$normbin[newdat$normalized >50 & newdat$normalized<=100]= "51-100"
newdat$normbin[newdat$normalized >100 & newdat$normalized<500]= "101-500"
newdat$normbin[newdat$normalized >500 & newdat$normalized<=1000]= "501-1000"
newdat$normbin[newdat$normalized >1000 & newdat$normalized<=1500]= "1001-1500"
newdat$normbin[newdat$normalized >1500]= ">1500"


# Add a column with labels so that we can just use these for Figures 2, S2-S3 when we create it in GIS
newdat$labels <- paste(newdat$prov, " (n=",newdat$numfish,")", sep = "")

# Write the data into a .csv file with all of the organized, binned variables.
#write.csv(newdat, "Longhurst_FishSummaryData_fullbinned_updated.csv") #(Map w/ all data; Fig S2A)

# if have run this 
#write.csv(newdat, "Longhurst_FishSummaryData_fullbinned_wMP.csv") #(Map w/ Microplastics only; Fig 2)
write.csv(newdat, "Longhurst_FishSummaryData_fullbinned_M3.csv") #(Map w/ Method 3 only; Fig S2B)

#### To bring this file into QGIS for mapping:
# 1. Import as a vector layer;
# 2. Join this layer by province code with the shape file for the Longhurst provinces that can be downloaded from here: http://www.marineregions.org/downloads.php#longhurst.
# 3. Adjust colors, labels, and variables using the "Properties" menu in QGIS.

# Figure 2 was made by color-coding the Longhurst province shape file according to the "aveplastbin" column from the microplastic only data. The labels were assigned from the "labels" column. 

#Figure S2 was made from a) the full dataset, and b) by subsetting the data to studies that used method 3. The labels were assigned from the "labels" column. 

#Figure S3 was made by color coding from the full data set according to the "numstudiesbin" column.


# Figure 3, plastic ingestion by depth and habitat ----
d_full %>% 
  pull(average_depth) %>% 
  quantile(c(0.025, 0.25, 0.5, 0.85, 0.90, 0.95, 0.99), na.rm = TRUE)


Fig_3 <- ggplot(
  filter(d_full, includes_microplastic == "Y", Found != "NA"), 
  aes(average_depth, prop_w_plastic, size = N)) +
  geom_point(aes(color = prop_w_plastic), alpha = 0.6) + 
  geom_smooth(aes(weight = N), col = "gray 30", method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Found, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_x_reverse() +
  scale_color_gradientn(colours = c("steelblue4",
                                    "darkgoldenrod1",
                                    "darkorange", "orangered1",
                                    "firebrick1", "red3", "red4"), 
                        name = "Proportion with \ningested plastic") +
  xlim(350,0) +
  xlab("Average depth (m)") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic(base_size = 16)
Fig_3 + guides(size = FALSE, color = FALSE)

dev.copy2pdf(file="Fig_AD_F.pdf", width=6, height=8)





# Figure 4 risk plot  ---- 
risk_plot <- d_sp_sum %>% 
  drop_na(commercial) %>% 
  ggplot(aes(log10(Sample_size), Sp_mean)) +
  geom_point(aes(color = Sp_mean, shape = commercial, size = studies_cat), 
             alpha = 0.8) +
  geom_hline(yintercept = 0.25, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = log10(10), linetype="dashed", color = "grey50") +
  xlab("Log[Sample Size]") +
  ylab("Species-specific plastic ingestion incidence (FO)") +
  labs(shape = "Commercial status", shape = "Commercial Status", size =  "Number of studies") +
  scale_color_gradientn(colours = c("steelblue4",
                                    "darkgoldenrod1",
                                    "darkorange", "orangered1",
                                    "firebrick1", "red3", "red4"), 
                        name = "Proportion with \ningested plastic") +
  scale_size_continuous(breaks = seq(from = 1, to = 3, by = 1), 
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n>3)"),
                        range = c(1.5, 5)) +
  annotate("text", x = c(0.4, 3.5, 0.4, 3.5),
           y=c(0.8, 0.8, 0.0725, 0.0725),
           label = c("high incidence, data poor", "high incidence, data rich",
                     "low incidence, data poor", "low incidence, data rich")) +
  theme_classic(base_size = 16) +
  guides(shape = guide_legend(override.aes = list(size = 3)))
risk_plot

dev.copy2pdf(file="risk_plot_SciAd.pdf", width=12, height=7)


# Figure 5, FO over time (2010-present), and species accumulation curves----
FO_year_2010 <- d_full %>% 
  drop_na(N,prop_w_plastic, publication_year) %>% 
  filter(publication_year >2009, method_type == 3) %>%   #; CAN TOGGLE THIS IN AND OUT; BUT NOW IT"S IN D FULL TOO
  mutate(commercial = fct_collapse(commercial,
                            Commercial = c("commercial", "highly commercial"),
                            Minor = c("minor commercial", "subsistence"),
                            None = "none")) 

FO_year_lm <- lm(prop_w_plastic~publication_year, weights = N, data = FO_year_2010)
summary(FO_year_lm)

FO_year <- ggplot(data = FO_year_2010, 
       aes(publication_year,prop_w_plastic, weight = N, size= N)) +
  geom_point(aes(color = prop_w_plastic, shape = commercial), alpha = 0.6) +
  geom_smooth(method = "lm", show.legend = FALSE, size = 0.75, color = "black") +
  xlim(2009,2020) + 
  geom_hline(yintercept = 0.26, linetype="dashed", color = "grey50") +
  scale_color_gradientn(colours = c("steelblue4",
                                    "darkgoldenrod1",
                                    "darkorange", "orangered1",
                                    "firebrick1", "red3", "red4")) +
  scale_size_continuous(breaks = c(1, 10, 100, 500, 1000)) +
  scale_x_continuous(breaks=seq(2010, 2020, 1)) +
  scale_shape_manual(na.translate = F, values = shapes) +
  labs(x = "Publication year",
       y = "Proportion with ingested plastic",
       size = "Sample size") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
FO_year

dev.copy2pdf(file="FO_by_year_Method3.pdf", width=6.75, height=6)


# average FO in 2019
FO_year_2019 <- FO_year_2010 %>% filter(publication_year == 2019)
weighted.mean(FO_year_2019$prop_w_plastic, w = FO_year_2019$N)



## Cumulative unique entries in list of vectors
cum_unique <- function(l) {
  result <- integer(length(l))
  for (i in 1:length(l)) {
    result[i] <- length(unique(unlist(l[1:i])))
  }
  result
}
d_rarefaction_all <-  d_full %>%
  filter(method_type == 3) %>%  # Can TOGGLE in and out
  group_by(publication_year) %>% 
  summarize(species = list(unique(binomial)),
            annual_N = sum(N, na.rm = TRUE)) %>% 
  ungroup %>% 
  mutate(cum_species = cum_unique(species),
         cum_n = cumsum(annual_N),
         cpue = cum_species / cum_n)
d_rarefaction_ingest_only <-  d_full %>%  
  filter(prop_w_plastic > 0) %>% 
  group_by(publication_year) %>% 
  summarize(species = list(unique(binomial)),
            annual_N = sum(N, na.rm = TRUE)) %>% 
  ungroup %>% 
  mutate(cum_species = cum_unique(species),
         cum_n = cumsum(annual_N),
         cpue = cum_species / cum_n)

rarefaction_plot <- ggplot() +
  geom_line(data = d_rarefaction_all, aes(cum_n, cum_species), color = "dark blue") +
  geom_line(data = d_rarefaction_ingest_only, aes(cum_n, cum_species), color = "dark red") +
  geom_point(data = d_rarefaction_all, aes(cum_n, cum_species)) +
  geom_point(data = d_rarefaction_ingest_only, aes(cum_n, cum_species)) +
  geom_text_repel(data = d_rarefaction_all, 
                  aes(cum_n, cum_species, label = publication_year), 
                  nudge_x = -1500, nudge_y = 10,
                  segment.color = "black") +
  geom_text_repel(data = d_rarefaction_ingest_only, 
                  aes(cum_n, cum_species, label = publication_year), 
                  nudge_x = 1500, nudge_y = -10,
                  segment.color = "black") +
  labs(x = "Cumulative number of individuals sampled",
       y = "Cumuluative number of species sampled") + 
  theme_classic(base_size = 16)
rarefaction_plot

dev.copy2pdf(file="Supp_rarefaction_plot.pdf", width=6, height=6)



# Supplemental figures----
# Figure S1, number of studies over time----
study_hist <- d %>% 
  group_by(publication_year, includes_microplastic) %>% 
  summarize(n_studies = n_distinct(source)) %>% 
  ggplot(aes(publication_year, n_studies)) + 
  geom_bar(aes(fill = includes_microplastic), stat = "identity") + 
  scale_fill_manual(labels = c("No", "Yes"), 
                    values = c("tomato3", "gray28")) +
  geom_smooth(se = FALSE) +
  theme_classic(base_size = 18) +
  labs(fill = "Recorded microplastics?",
       x = "Publication year",
       y = "Number of studies") +
theme(axis.title.x = element_blank())
study_hist 


study_hist_MT <- d %>% 
  group_by(publication_year, method_type) %>% 
  summarize(n_studies = n_distinct(source)) %>% 
  mutate(method_type = factor(method_type)) %>% 
  ggplot(aes(publication_year, n_studies)) + 
  geom_bar(aes(fill = method_type), stat = "identity") + 
  scale_fill_manual(labels = c("1", "2", "3"), 
                    values = c("darkslateblue",  "darkgoldenrod", "darkslategray")) +
  geom_smooth(se = FALSE, color = "gray48") +
  theme_classic(base_size = 18) +
  labs(fill = "Method",
       x = "Publication year",
       y = "Number of studies") +
  theme(axis.title.x = element_blank())
study_hist_MT 

dev.copy2pdf(file="studies_by_year_MT.pdf", width=12, height=7)

study_hist_AW <- d %>% 
  group_by(publication_year, adjacency_water) %>% 
  summarize(n_studies = n_distinct(source)) %>% 
  ggplot(aes(publication_year, n_studies)) + 
  geom_bar(aes(fill = adjacency_water), stat = "identity") + 
  scale_fill_manual(values = c( "burlywood4", "seagreen4", "dodgerblue4")) +
  geom_smooth(se = FALSE, color = "gray25") +
  theme_classic(base_size = 18) +
  labs(fill = "Region",
       x = "Publication year",
       y = "Number of studies") 
study_hist_AW 

dev.copy2pdf(file="studies_by_year_AW.pdf", width=12, height=7)

Study_hist_comb <- ggarrange(study_hist, study_hist_MT, study_hist_AW, 
                                labels = c("A","B", "C"), 
                             heights = c(1,1,1.1),
                                font.label = list(size = 20),
                                legend = "none",
                                ncol = 1, nrow = 3)
Study_hist_comb

dev.copy2pdf(file="Study_hist_comb.pdf", width=6, height=10)



# Figure S2.---- 
  ##Reference code for Figure 2.

# Figure S3.---- 
  ##Reference code for Figure 2.

# Figure S4. Mean plastic concentration per province ----
# load packages and data
library(sf)
library(rgdal)
library(raster)
library(maptools)

#note that to read in the shape file you need the full address (https://datacarpentry.org/r-raster-vector-geospatial/06-vector-open-shapefile-in-r/)

lh_prov <- st_read("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/longhurst_v4_2010/Longhurst_world_v4_2010.shp") # read in the Longhurst province shape file (referenced in Figures 2 and S4)
crs(lh_prov)
lh_prov.2 <- fortify(lh_prov) #fortify makes this into a data frame object for ggplot2

# We use this file to first assign adjacency values that will ultimately be used in the model: 1 for LH province touching major landmass (not islands), 0 for not touching major landmass

lh_prov.2$adjacency <- ifelse(lh_prov.2$ProvCode=="ALSK"| lh_prov.2$ProvCode=="CCAL"|lh_prov.2$ProvCode=="CAMR"|lh_prov.2$ProvCode=="BERS"|lh_prov.2$ProvCode=="BPLR"|lh_prov.2$ProvCode=="CARB"|lh_prov.2$ProvCode=="CHIL"|lh_prov.2$ProvCode=="FKLD"|lh_prov.2$ProvCode=="BRAZ"|lh_prov.2$ProvCode=="APLR"|lh_prov.2$ProvCode=="EAFR"|lh_prov.2$ProvCode=="BENG"|lh_prov.2$ProvCode=="GUIN"|lh_prov.2$ProvCode=="NWCS"|lh_prov.2$ProvCode== "CNRY"|lh_prov.2$ProvCode=="MEDI"|lh_prov.2$ProvCode== "NECS"|lh_prov.2$ProvCode=="SARC"|lh_prov.2$ProvCode=="REDS"|lh_prov.2$ProvCode== "ARAB"|lh_prov.2$ProvCode== "INDW"|lh_prov.2$ProvCode== "INDE"|lh_prov.2$ProvCode== "AUSW"|lh_prov.2$ProvCode== "BERS"|lh_prov.2$ProvCode=="SUND"|lh_prov.2$ProvCode== "AUSE"|lh_prov.2$ProvCode== "KURO"|lh_prov.2$ProvCode== "CHIN"|lh_prov.2$ProvCode=="TASM"|lh_prov.2$ProvCode=="NEWZ"|lh_prov.2$ProvCode=="SPSG"|lh_prov.2$ProvCode== "ARCH", 1, 0)
sum(lh_prov.2$adjacency)

# to calculate average plastic concentration per region, we need to bring in open-access data from Van Sebille et al. 2015. Here we use the Van Sebille model of abundance, presented in the format of three .csvs: one with abundance, one with latitude, and one with longitude. 
abund <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/vansebillemodel_abundance.csv")
poll_lat <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/latitudes.csv", header = FALSE)
poll_lon <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/longitudes.csv", stringsAsFactors = FALSE)

#Creating new object for plastic rasters
# clean up longitude values
x_vals <- as.numeric(gsub("X", "", colnames(poll_lon))) #gsub = regular expression that will drop the X from in front of the numeric values in our x data
x_vals[x_vals > 179]  = x_vals[x_vals > 179] - 360 #need to shift the view on the map, so essentially flip the map around the 180.
# incorporating into raster object
raster_obj <- list(z = as.matrix(abund)[, order(x_vals)], 
                   x = x_vals, 
                   y = as.numeric(poll_lat$V1))

# Create a new raster() file
poll_raster <- raster(x = raster_obj$z, # Matrix values of plastic
                      
                      # Defining endpoints of the raster
                      xmn = min(raster_obj$x),
                      xmx = max(raster_obj$x),
                      ymn = min(raster_obj$y),
                      ymx = max(raster_obj$y),
                      
                      # Setting coordinate reference system to be equivalent to longhurst
                      crs = crs(lh_prov))
# checking plot function
plot(st_geometry(lh_prov), add = TRUE, fill = NULL)
### Next step: bring out average attribute per polygon
extracted_vals <- extract(poll_raster, lh_prov) #this extracts values in the poll_raster per polygon
str(extracted_vals) #this should give us a list of all of the values in each of 54 different LH objects
mean_poll_abund<- unlist(lapply(extracted_vals, mean, na.rm=T)) #this should give us the mean values for each of the 54 LH objects
lapply(extracted_vals, range, na.rm=T) #note that these units are #/km^2

fullmapdat<- cbind(lh_prov.2, mean_poll_abund) #pairing each province with the appropriate mean pollution concentration values 

### need to bin these vals for the map - values are AVERAGE #/km^2
summary(mean_poll_abund)
fullmapdat$plasticbin <- rep(NA, length(fullmapdat$mean_poll_abund))
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=1 & fullmapdat$mean_poll_abund <10]= "10e0"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=10 & fullmapdat$mean_poll_abund <100]= "10e1"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=100 & fullmapdat$mean_poll_abund <1000]= "10e2"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=1000 & fullmapdat$mean_poll_abund <10000]= "10e3"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=10000 & fullmapdat$mean_poll_abund <100000]= "10e4"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=100000 & fullmapdat$mean_poll_abund <1000000]= "10e5"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=100000 & fullmapdat$mean_poll_abund <1000000]= "10e6"
fullmapdat$plasticbin[fullmapdat$mean_poll_abund >=1000000]= ">10e6"
fullmapdat$plasticbin

# now we make make sure this a dataframe, rather than a spatial object
fullmap.df <- fortify(fullmapdat)
fullmap.df <- st_drop_geometry(fullmap.df) #need to remove geometry in order to get an actual dataframe (zero spatial component)
# to save this file as a .csv
write.csv(fullmap.df, "Spatial Information_microplastics.csv")

### This .csv file can then be imported into QGIS as a vector layer to generate color-coded maps (see instructions in code for Figure 2).

# Figure S5----
d_full_wo_gaps_AD_F <- d %>%
  filter(includes_microplastic == "Y") %>% 
  #filter(method_type == 3) %>%  # Can TOGGLE in and out
  drop_na(average_depth, Found, NwP, N)

glmm_FwP_AD_F <- glmer(cbind(NwP, N-NwP) ~ 
                         scale(average_depth)*Found + 
                         (1|order) + (1|source) + (1|method_type),  
                       na.action = "na.fail",
                       data = d_full_wo_gaps_AD_F, 
                       family = binomial)
summary(glmm_FwP_AD_F)
r.squaredGLMM(glmm_FwP_AD_F)


plot_model(glmm_FwP_AD_F, type = "re",
           grid = FALSE, 
           sort.est = "sort.all")[[2]] +  # Change between 1,2,3 to produce plots for each random effect
  ggtitle("") +
  #xlab("Order") +
  ylab("Random intercept") +
  #ylim(0.1,10) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1)

dev.copy2pdf(file="Random effect of order_GLMM.pdf", width=7, height=12)






# Possible supplemental figure: Species-level phylogeny ----

# building the basic tree 

breaks <- c(seq(1,nrow(d_sp_sum),50),nrow(d_sp_sum)+1)  # why are we doing this?

for (i in 1:(length(breaks)-1)){
  taxa <- as.character(d_sp_sum$binomial[breaks[i]:(breaks[i+1]-1)])
  taxa <- taxa[taxa != "" & !is.na(taxa)]
  
  resolved_namest <- tnrs_match_names(taxa)                          # I think this is where all the extra species fall out
  resolved_namest <- resolved_namest[!is.na(resolved_namest$unique_name),]
  if (i==1){
    resolved_namess <- resolved_namest
  } else {
    resolved_namess <- rbind(resolved_namess, resolved_namest)
  }
}
resolved_names <- resolved_namess
resolved_names <- resolved_names[resolved_names$flags!="INCERTAE_SEDIS_INHERITED",]

#plots species
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id, label_format = "name")

my_tree$tip.label<-gsub("_"," ",my_tree$tip.label) # removes underscore between genus and species names
my_tree$tip.label<-str_extract(my_tree$tip.label, "[A-Z][a-z]+ [a-z]+")

my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

#View(my_tree)

#my_tree2 = phylo4d(my_tree, d_sp_sum)


# first plot try, fan layout
p <- ggtree(my_tree, layout="fan", open.angle=0) + 
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  ggplot2::xlim(-0.6, 1.3) 
p

my_tree$tip.label <- as.factor(my_tree$tip.label)

# SUPPLEMENTAL PHYLOGENY FIGURE, ALL SPECIES
p %<+% d_sp_sum + 
  aes(color = Sp_mean) +
  geom_tiplab2(aes(label = paste0("italic('", label, "')"),
                   color=Sp_mean, angle = angle), parse = TRUE,
               size = 6, align = FALSE, hjust = -0.05) +
  geom_tippoint(aes(color = Sp_mean, shape = commercial, size = studies_cat)) +
  scale_color_gradientn(colours = c("steelblue4",
                                    "darkgoldenrod1",
                                    "darkorange", "orangered1",
                                    "firebrick1", "red3", "red4"),
                        name = "Proportion with \ningested plastic") +
  #scale_size(range = c(3, 7)) +
  scale_size_continuous(guide = FALSE, range = c(3, 7)) +
  scale_shape_discrete(na.translate = F) +
  labs(shape = "Commercial \nstatus") +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(3.5, "cm"),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.box = "horizontal") +
  guides(
    size = FALSE, shape = FALSE, color  = FALSE
    #shape = guide_legend(override.aes = list(size = 5))
  )


# save plots
dev.copy2pdf(file="Fish_plastic_phylo_d_mp_subset.pdf", width=50, height=50)
ggsave("Prelim_phylo_ggtree2.tiff", width = 42, height = 42, units = "in")
ggsave("Prelim_phylo_ggtree2.eps", width = 45, height = 45, units = "in")
ggsave("Fish_plastic_phylo_d_mp_subset_nolegend.pdf", width = 45, height = 45, units = "in")





