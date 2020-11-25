library(tidyverse)
library(sf)
library(tidyverse)
library(gbm)
library(dismo)
library(mgcv)
library(lme4)
library(gamm4)
library(readxl)
library(readr)
library(ggplot2)
library(rnaturalearth)
library(rgeos)
library(RColorBrewer)
library(raster)
library(dplyr)
getwd()
#longhurst <- read_sf("/Users/agmc/Box/Microplastics/Fish-plastics-meta-analysis-code_AGM/Complete_files_Savoca et al/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
longhurst <- read_sf("Longhurst_world_v4_2010.shp")

#### Add Microplastic data - taken from previous code####
data <- read.csv("Plastics ingestion records fish master_final_GCB_v2.csv")
# bring in reliability scores from New_test_code
data <- data %>% filter(Includes.microplastic == "Y"|Includes.microplastic == "Y?") # include this line if you're looking for studies that have to include microplastics
nrow(data)
# NEED TO REMOVE ESTUARINE STUDIES ##
data <- data %>% filter(Water.type == "marine")

# Additional edits here for R1
data <- data %>% filter(Source != "Sun et al. 2019")  #duplicate with Zhao et al. 2019

data <- data[!is.na(data$Oceanographic.province..from.Longhurst.2007.) &data$Oceanographic.province..from.Longhurst.2007.!="",]
nrow(data)
# want to get average proportion of plastic per province
summary(data)
head(data)
# data of interest
data2 <- data[,c("Binomial", "Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "NwP", "N", "Source")]
head(data2)
colnames(data2) <- c("Species", "OceanProv", "PropPlastic", "NwP", "N", "Source")
head(data2)
length(table(data2$OceanProv))
data3<- data2[order(data2$OceanProv),]
head(data3)

data3$OceanProv <- as.factor(data3$OceanProv)
levels(data3$OceanProv) <- c(levels(data3$OceanProv), "CHIL", "BPLR", "NASE", "NPPF") # need to change these to match the publicly available shape file
data3$OceanProv[data3$OceanProv=="BPRL"] <- "BPLR"
data3$OceanProv[data3$OceanProv=="HUMB"] <- "CHIL"
data3$OceanProv[data3$OceanProv=="NAST E"] <- "NASE"
data3$OceanProv[data3$OceanProv=="NPSE"] <- "NPPF"

table(data3$OceanProv)
data3$OceanProv <- droplevels(data3$OceanProv)
levels(data3$OceanProv)
#######################################################################################
#################### IF ALREADY HAVE FULL CLEANED DATA (d_full_R1) #################### 
data3 <- as.data.frame(d_full_R1) #fom New_test_code

### set up new dataframe 
prov <- unique(data3$ProvCode)
prov
length(prov)
aveplast <- rep(NA, length(prov))
numfish <- rep(NA, length(prov))
numstudies <- rep(NA, length(prov))
numspecies <- rep(NA, length(prov))
normalized <- rep(NA, length(prov))
medianreliability <- rep(NA, length(prov))

 ## using this to double check the produced averages
sub <- data3[data3$ProvCode=="KURO",]
sub
unique(sub$source)
length(unique(sub$source))
overallprop <- sum(sub$NwP, na.rm=T)/sum(sub$N, na.rm=T)
overallprop
normalized <- overallprop*sum(sub$N, na.rm=T)
normalized
sum(sub$N, na.rm=T)
mean(sub$prop_w_plastic, na.rm=T) #3 = .1077803
sum(sub$N)
nrow(sub)
length(unique(sub$species))
# gives us reliability for each study
rels <- sub %>% group_by(source) %>% summarise(reliability=mean(overall_reliability))
median(rels$reliability)

for (i in 1:length(prov)){
  sub <- data3[data3$ProvCode==prov[i],]
  aveplast[i] <- sum(sub$NwP, na.rm=T)/sum(sub$N, na.rm=T)
  numfish[i] <- sum(sub$N, na.rm=T)
  numstudies[i] <- length(unique(sub$source))
  numspecies[i] <- length(unique(sub$species))
  normalized[i] <- aveplast[i]*sum(sub$N, na.rm=T)
  rels <- sub %>% group_by(source) %>% summarise(reliability=mean(overall_reliability))
  medianreliability[i] <- median(rels$reliability)
}

aveplast
length(aveplast)
numfish
table(numfish)
numstudies
head(data3)
numspecies
medianreliability

newdat <-  data.frame(prov, aveplast, numfish, numstudies, numspecies, normalized, medianreliability)
newdat


## binning data
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

head(newdat$aveplastbin, 50)
sort(newdat$aveplastbin)
# number of studies
summary(newdat$numstudies)
newdat$numstudiesbin <- rep(NA, length(prov))
newdat$numstudiesbin[newdat$numstudies == 1] = "1"
newdat$numstudiesbin[newdat$numstudies >1 & newdat$numstudies<=5]= "2-5"
newdat$numstudiesbin[newdat$numstudies >5 & newdat$numstudies<=10]= "6-10"
newdat$numstudiesbin[newdat$numstudies >=10]= ">10"

head(newdat, 50)

# number of fish/study
summary(newdat$numfish)
newdat$numfishbin <- rep(NA, length(prov))
newdat$numfishbin[newdat$numfish >1 & newdat$numfish<10]= "< 10"
newdat$numfishbin[newdat$numfish >=10 & newdat$numfish<=50]= "10-50"
newdat$numfishbin[newdat$numfish >50 & newdat$numfish<=100]= "51-100"
newdat$numfishbin[newdat$numfish >100 & newdat$numfish<500]= "101-500"
newdat$numfishbin[newdat$numfish >500 & newdat$numfish<=1000]= "501-1000"
newdat$numfishbin[newdat$numfish >1000 & newdat$numfish<=1500]= "1001-1500"
newdat$numfishbin[newdat$numfish >1500]= ">1500"

head(newdat, 50)

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

head(newdat, 50)
newdat$aveplast
sort(newdat$numstudies)
sort(newdat$numfish)
head(newdat)

newdat$labels <- paste(newdat$prov, " (n=",newdat$numfish,")", sep = "")

head(newdat)

#write.csv(newdat, "Longhurst_FishSummaryData_fullbinned_wMP.csv")
#write.csv(newdat, "Longhurst_FishSummaryData_fullbinned.csv")

############################ Spatial Data Analysis/Mapping #################################
library(sf)
library(rgdal)
library(raster)
library(maptools)
#note that to read in the shape file you need the full address (https://datacarpentry.org/r-raster-vector-geospatial/06-vector-open-shapefile-in-r/)
# ggplot2 will only work with a data.frame object

lh_prov <- longhurst
lh_prov
crs(lh_prov)
lh_prov.2 <- fortify(lh_prov) #fortify makes this into a data frame object

## adding adjacency values: 1 for LH province touching major landmass (not islands), 0 for not touching major landmass

lh_prov.2$adjacency <- ifelse(lh_prov.2$ProvCode=="ALSK"| lh_prov.2$ProvCode=="CCAL"|lh_prov.2$ProvCode=="CAMR"|lh_prov.2$ProvCode=="BERS"|lh_prov.2$ProvCode=="BPLR"|lh_prov.2$ProvCode=="CARB"|lh_prov.2$ProvCode=="CHIL"|lh_prov.2$ProvCode=="FKLD"|lh_prov.2$ProvCode=="BRAZ"|lh_prov.2$ProvCode=="APLR"|lh_prov.2$ProvCode=="EAFR"|lh_prov.2$ProvCode=="BENG"|lh_prov.2$ProvCode=="GUIN"|lh_prov.2$ProvCode=="NWCS"|lh_prov.2$ProvCode== "CNRY"|lh_prov.2$ProvCode=="MEDI"|lh_prov.2$ProvCode== "NECS"|lh_prov.2$ProvCode=="SARC"|lh_prov.2$ProvCode=="REDS"|lh_prov.2$ProvCode== "ARAB"|lh_prov.2$ProvCode== "INDW"|lh_prov.2$ProvCode== "INDE"|lh_prov.2$ProvCode== "AUSW"|lh_prov.2$ProvCode== "BERS"|lh_prov.2$ProvCode=="SUND"|lh_prov.2$ProvCode== "AUSE"|lh_prov.2$ProvCode== "KURO"|lh_prov.2$ProvCode== "CHIN"|lh_prov.2$ProvCode=="TASM"|lh_prov.2$ProvCode=="NEWZ"|lh_prov.2$ProvCode=="SPSG"|lh_prov.2$ProvCode== "ARCH", 1, 0)
sum(lh_prov.2$adjacency)

lh_prov.2$adjacency
head(lh_prov.2)


######### Base map #### 
# World map
world_map <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))
# Base map
kk <- ggplot() +
  geom_sf(data = world_map, size = .2, fill = "gray30", col = "gray90") +
  borders(colour="NA", fill="darkgray")  +
  theme(panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(color = "white", linetype = "dashed", size = 0.5))
kk

# simplify the object to make it 'usable'
longhurst <- lh_prov %>% 
  sf::st_simplify(dTolerance = 0.01) %>% 
  dplyr::group_by(ProvCode,ProvDescr) %>% 
  dplyr::summarise()
plot(longhurst)

# plot
kk+  
  geom_sf(data = lh_prov, aes(fill = ProvCode), size = .2, col = "grey50", alpha=.4)+
  ggtitle(paste("Longhurst Biogeochemical Provinces -", length(unique(longhurst$ProvCode)),"provinces"))+
  theme(legend.position="none")+
  geom_sf_text(data = longhurst %>% group_by(ProvCode) %>% dplyr::summarise(n()), aes(label = ProvCode), colour = "grey20", check_overlap=TRUE)+
  coord_sf(expand = FALSE)


######### Average Plastic Pollution Map ################
abund <- read.csv("vansebillemodel_abundance.csv")
poll_lat <- read.csv("latitudes.csv", header = FALSE)
poll_lon <- read.csv("longitudes.csv", stringsAsFactors = FALSE)

####### Creating new object for plastic rasters
head(abund)
head(poll_lat)
head(poll_lon)

x_vals <- as.numeric(gsub("X", "", colnames(poll_lon))) #gsub = regular expression that will drop the X from in front of the numeric values in our x data (also, lon reads in as column names)
x_vals[x_vals > 179]  = x_vals[x_vals > 179] - 360 #need to shift the view on the map, so essentially flip the map around the 180. Probably need to get rid of 360 because 0 and 360 are redundant


###### IF YOU WANT TO GET RID OF THE LINE AT 0/360 or 0/0 (if we are looking at 0-180 scale), run this code instead ###
# Creating raster list to generate full object 
raster_obj <- list(z = as.matrix(abund)[, order(x_vals[-1])], #list with attributes as z, and coordinates as x and y
    x = x_vals[-1], #got rid of the initial 0 value, because 0 and 360 (or in this case, if we are going between -180 and 180, 0 and 0) overlap. This gets rid of the weird disconnect in the middle of the map
    y = as.numeric(poll_lat$V1))
##############################################################################################################

#to calculate average plastic, will keep the weird line in because we want the plastic values at 0 and 360
raster_obj <- list(z = as.matrix(abund)[, order(x_vals)], #list with attributes as z, and coordinates as x and y
                   x = x_vals, #got rid of the initial 0 value, because 0 and 360 (or in this case, if we are going between -180 and 180, 0 and 0) overlap. This gets rid of the weird disconnect in the middle of the map
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

# Plotting Z on the log10 scale (Van Sebille et al (2015, ERL) paper); see that this works; don't want calculations on this scale
plot(poll_raster)

plot(st_geometry(lh_prov), add = TRUE, fill = NULL)

### Next step: bring out average attribute per polygon
extracted_vals <- raster::extract(poll_raster, lh_prov) #this extracts values in the poll_raster per polygon
str(extracted_vals) #this should give us a list of all of the values in each of 54 different LH objects
mean_poll_abund<- unlist(lapply(extracted_vals, mean, na.rm=T)) #this should give us the mean values for each of the 54 LH objects
lapply(extracted_vals, range, na.rm=T) #note that these units are #/km^2

fullmapdat<- cbind(lh_prov.2, mean_poll_abund)
fullmapdat

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

head(fullmapdat) #full dataset with all information, including geometries

# to save this information in dataframe form
fullmap.df <- fortify(fullmapdat)
fullmap.df <- st_drop_geometry(fullmap.df) #need to remove geometry in order to get an actual dataframe (zero spatial component)
head(fullmap.df)
nrow(fullmap.df)

full_sampled_data <- merge(fullmap.df, newdat, by.x="ProvCode", by.y="prov", sort = TRUE)

#write.csv(fullmap.df, "Spatial Information_microplastics.csv")
write.csv(full_sampled_data, "Spatial Information_onlysampled.csv")

# to plot this data on map
colors = brewer.pal(6,"Reds")

pollution_map <- kk+  
  geom_sf(data = fullmapdat, aes(fill = log10(mean_poll_abund)), size = .2, alpha=.7)+
  scale_fill_gradientn(colours=colors) + 
  labs(x="", y="") + 
  ggtitle(paste("Mean Pollution per LH Prov -", length(unique(fullmapdat$ProvCode)),"provinces"))+
  ggrepel::geom_label_repel(data=longhurst %>% group_by(ProvCode) %>% dplyr::summarize(n()), aes(label = ProvCode, geometry=geometry, fontface="bold"), label.padding = unit(0.1, "lines"), size=2.5,
                            stat = "sf_coordinates",
                            min.segment.length = 0,
                            label.size = NA, 
                            color = "black") +
  theme(axis.title=element_text(size=14, face="bold"), plot.title = element_text(hjust=.5, size=16, face="bold"), legend.title=element_blank())

pollution_map
#check and see if this makes sense
fullmap.df[order(fullmap.df$mean_poll_abund),]


###### Maps based on our data ######
head(newdat)

lh_prov.2
data_map <- merge(lh_prov.2, newdat, by.x="ProvCode", by.y="prov", all.x=T)
#data_map <- merge(fullmapdat, newdat, by.x="ProvCode", by.y="prov", all.x=T)
data_map

colorvals = c("steelblue4",
                        "darkgoldenrod1",
                        "darkorange", "orangered1",
                        "firebrick1", "red3", "red4", "white")

data_map$aveplastbin <- as.factor(data_map$aveplastbin)
data_map$aveplastbin <- ordered(data_map$aveplastbin, levels = c("NA", "<.10", ".11-.20", ".21-.30", ".31-.40", ".41-.50", ".71-.80", "1"))

tail(data_map)
###### map for average plastic binned (MP map only)

# map for average plastic binned (MP map only)
plastic_map2 <- kk +
  
  geom_sf(data = data_map,aes(fill = aveplastbin), size = .2, color="grey", alpha=.8) +
  #borders(colour=NA, fill = "white") +
  scale_fill_manual(values=colorvals) +
  #scale_fill_manual(values=colors(length(unique(data_map$aveplastbin))))
  labs(x="", y="") +
  ggrepel::geom_label_repel(data=data_map %>% 
                              group_by(labels) %>% 
                              dplyr::summarize(n()), 
                            aes(label = labels, geometry=geometry, fontface="bold"), 
                            label.padding = unit(0.1, "lines"), size=2.5,
                            stat = "sf_coordinates",
                            min.segment.length = 0,
                            label.size = NA, 
                            color = "black") +
  theme(panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(color = "white", 
                                        linetype = "dashed", size = 0.5),
        axis.title=element_text(size=14, face="bold"), 
        plot.title = element_text(hjust=.5, size=16, face="bold"), 
        legend.title=element_blank())

plastic_map2

dev.copy2pdf(file="Fig_2_allMP.pdf", width=12, height=8)

# map for number of studies
colorvals2 <-  brewer.pal(4,"YlGn")

data_map$numstudiesbin <- as.factor(data_map$numstudiesbin)
data_map$numstudiesbin <- ordered(data_map$numstudiesbin, levels = c("NA", "1", "2-5", "6-10", ">10"))

#data_map$labels <- ifelse(!is.na(data_map$labels), data_map$labels, "")

studies_map <- kk+  
  geom_sf(data = data_map, aes(fill = numstudiesbin), color="grey", size = .2, alpha=.6)+
  scale_fill_manual(values=colorvals2) +
  #scale_fill_manual(values=colors(length(unique(data_map$aveplastbin))))
  ggtitle("Number of studies per region")+  
  labs(x="", y="") +
  ggrepel::geom_label_repel(data=data_map %>% group_by(labels) %>% dplyr::summarize(n()), aes(label = labels, geometry=geometry, fontface="bold"), label.padding = unit(0.1, "lines"), size=2.5,
                            stat = "sf_coordinates",
                            min.segment.length = 0,
                            label.size = NA, 
                            color = "black") +
  theme(axis.title=element_text(size=14, face="bold"), plot.title = element_text(hjust=.5, size=16, face="bold"), legend.title=element_blank())


studies_map

## map for number of fish sampled
colorvals3 <-  brewer.pal(7,"RdPu")

data_map$numfishbin <- as.factor(data_map$numfishbin)
data_map$numfishbin <- ordered(data_map$numfishbin, levels = c("NA", "< 10", "10-50", "51-100", "101-500", "501-1000","1001-1500", ">1500"))


### with code from https://www.johan-rosa.com/2019/11/13/creating-maps-with-sf-and-ggplot2/
fish_map.2 <- kk+  
  geom_sf(data = data_map, aes(fill = numfishbin), color="grey", size = .2, alpha=.6)+
  scale_fill_manual(values=colorvals3) +
  #scale_fill_manual(values=colors(length(unique(data_map$aveplastbin))))
  ggtitle("Number of fish sampled per region")+  
  labs(x="", y="") +
  ggrepel::geom_label_repel(data=data_map %>% group_by(labels) %>% dplyr::summarize(n()), aes(label = labels, geometry=geometry, fontface="bold"), label.padding = unit(0.1, "lines"), size=2.5,
    stat = "sf_coordinates",
    min.segment.length = 0,
    label.size = NA, 
    color = "black") +
  theme(axis.title=element_text(size=14, face="bold"), plot.title = element_text(hjust=.5, size=16, face="bold"), legend.title=element_blank())

fish_map.2 

## map for reliability
colorvals4 <-  brewer.pal(8,"Blues")
data_map$medianreliability <- as.factor(data_map$medianreliability)
#data_map$medianreliability <- ordered(data_map$medianreliability, levels = c("NA", "< 10", "10-50", "51-100", "101-500", "501-1000","1001-1500", ">1500"))


reliability <- kk+  
  geom_sf(data = data_map, aes(fill = medianreliability), color="grey", size = .2, alpha=.6)+
  scale_fill_manual(values=colorvals4) +
  #scale_fill_manual(values=colors(length(unique(data_map$aveplastbin))))
  ggtitle("Median reliability score per region")+  
  labs(x="", y="") +
  ggrepel::geom_label_repel(data=data_map %>% group_by(labels) %>% dplyr::summarize(n()), aes(label = labels, geometry=geometry, fontface="bold"), label.padding = unit(0.1, "lines"), size=2.5,
                            stat = "sf_coordinates",
                            min.segment.length = 0,
                            label.size = NA, 
                            color = "black") +
  theme(axis.title=element_text(size=14, face="bold"), plot.title = element_text(hjust=.5, size=16, face="bold"), legend.title=element_blank())

reliability 
