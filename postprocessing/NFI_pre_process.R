#### Script to pre-process NFI data to fit format of LPJ-GUESS output.


library(dplyr)
library(DGVMTools)
setwd("C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\updated_data")
input_data <- read.table('tmt_grid_4evaluation_03-01-24.txt',header=TRUE)
colnames(input_data)

input_data[is.na(input_data)] <- 0


################## Prepare data to match DGVMTools GUESS format #########################
Lat_lon_year<-input_data[c("Lat","Lon","census.date", "country")] # add "country" to de-select or select certain countries later
Lat_lon_year$census.date<-as.integer(Lat_lon_year$census.date)
colnames(Lat_lon_year)[3] <- "Year"

# Optional filter by country
Filtered_NFI <- input_data
mask <- Filtered_NFI$country %in% c("Spain", "Germany", "Netherlands")
Filtered_NFI <- Filtered_NFI[mask, ]
write.table(Filtered_NFI,"C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\ref\\Filtered_NFI.txt")

## Create NFI .out for NPP
Forest_sum<-input_data[ , "npp2_mean"]
pann_mean<-cbind(Lat_lon_year,Forest_sum)
#mask <- pann_mean1$country %in% c("Germany") #commented lines are to deselect sertain countries
#pann_mean <- pann_mean1[!mask, ]
#pann_mean <- pann_mean[, !(names(pann_mean) == 'country')]
pann_mean[, 4]<-((pann_mean[,4]/10)/2)/0.8 #AGB conversion
write.table(pann_mean,"C:\\Users\\Admin\\Documents\\tellus\\Data_240122\\european_applications\\cmass_wood_inc_plus10cm_sts.out", sep = "\t",row.names = FALSE, quote = FALSE)

## Create NFI .out for Conifer Share
Total<-(input_data[ , 113])*100
biomass_conifer_share<-cbind(Lat_lon_year,Total)
write.table(biomass_conifer_share, "C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\ref\\biomass_conifer_share.out", sep = "\t",row.names = FALSE, quote = FALSE)

## Create NFI .out for Stem density
mean_stemdensity_bins<-input_data[, 12:36]
Forest_sum<-input_data[,"stemdensity_mean"]
dens_plus10cm_sts_ext_patch_mean1<-cbind(Lat_lon_year,mean_stemdensity_bins)
dens_plus10cm_sts_ext_patch_mean <- cbind(dens_plus10cm_sts_ext_patch_mean1, Forest_sum)
#mask <- dens_plus10cm_sts_ext_patch_mean2$country %in% c("Norway") #commented lines are to deselect sertain countries
#dens_plus10cm_sts_ext_patch_mean <- dens_plus10cm_sts_ext_patch_mean2[mask, ]
#dens_plus10cm_sts_ext_patch_mean <- dens_plus10cm_sts_ext_patch_mean[, !(names(dens_plus10cm_sts_ext_patch_mean) == 'country')]
#dens_plus10cm_sts_ext_patch_mean[, 4:30]<-dens_plus10cm_sts_ext_patch_mean[,4:30]
write.table(dens_plus10cm_sts_ext_patch_mean, "C:\\Users\\Admin\\Documents\\tellus\\Data_240122\\european_applications\\dens_plus10cm_sts.out", sep = "\t",row.names = FALSE, quote = FALSE)

# Test plot
ggplot(dens_plus10cm_sts_ext_patch_mean, aes(x=Lon, y=Lat, color=Forest_sum)) + 
  geom_point(size=2)

plotSpatial(dens_plus10cm_sts_ext_patch_mean, layers = "Forest_sum")

## Create NFI .out for Biomass
mean_biomass_ha<-input_data[ , 37:61]
Forest_sum<-input_data[,4]
cmass_plus10cm_sts_ext_patch_mean1<-cbind(Lat_lon_year,mean_biomass_ha)
cmass_plus10cm_sts_ext_patch_mean <- cbind(cmass_plus10cm_sts_ext_patch_mean1, Forest_sum)
mask <- cmass_plus10cm_sts_ext_patch_mean$country %in% c("Germany", "Spain","Netherlands") # Optional mask
cmass_plus10cm_sts_ext_patch_mean <- cmass_plus10cm_sts_ext_patch_mean[mask, ]
cmass_plus10cm_sts_ext_patch_mean <- cmass_plus10cm_sts_ext_patch_mean[, !(names(cmass_plus10cm_sts_ext_patch_mean) == 'country')]
cmass_plus10cm_sts_ext_patch_mean[, 4:29]<-((cmass_plus10cm_sts_ext_patch_mean[,4:29]/10)/2)/0.8 #1.6 AGB
write.table(cmass_plus10cm_sts_ext_patch_mean, "C:\\Users\\Admin\\Documents\\tellus\\Data_240122\\european_applications\\filtered_cmass.out", sep = "\t", row.names = FALSE, quote = FALSE)

## Create NFI .out for DBH
dbh_mean<-input_data[ , 4]
#Total<-rowSums(dbh_mean[,1:25])
diam_plus10cm_sts_ext_patch_mean<-cbind(Lat_lon_year,dbh_mean)
#diam_plus10cm_sts_ext_patch_mean <- cbind(diam_plus10cm_sts_ext_patch_mean1, Total)
write.table(diam_plus10cm_sts_ext_patch_mean, "E:\\Tellus\\ForMMI_forest_evaluation\\ref\\diam_plus10cm_sts_ext_patch_mean.out", sep = "\t", row.names = FALSE, quote = FALSE)

## Create NFI .out for QDBH
Forest_sum<-input_data[ , 10]
Forest_sum<-Forest_sum/10
diam_g_plus10cm_sts_ext_patch_mean<-cbind(Lat_lon_year,Forest_sum)
#diam_g_plus10cm_sts_ext_patch_mean <- cbind(diam_g_plus10cm_sts_ext_patch_mean1, Total)
write.table(diam_g_plus10cm_sts_ext_patch_mean, "C:\\Users\\Admin\\Documents\\tellus\\Data_240122\\european_applications\\diam_g_plus10cm_sts.out", sep = "\t", row.names = FALSE, quote = FALSE)

## Create NFI .out for tree total obs
obs_tree_total<-input_data[ , 87:111]
Total<-rowSums(obs_tree_total[,1:25])
diamstruct_forest1<-cbind(Lat_lon_year,obs_tree_total)
diamstruct_forest <- cbind(diamstruct_forest1, Total)
write.table(diamstruct_forest, "E:\\Tellus\\ForMMI_forest_evaluation\\ref\\diamstruct_forest.out", sep = "\t", row.names = FALSE, quote = FALSE)

## Create NFI .out for tree mean obs
obs_tree_mean<-input_data[ , 62:86]
Total<-rowSums(obs_tree_mean[,1:25])
diamstruct_forest1<-cbind(Lat_lon_year,obs_tree_mean)
diamstruct_forest <- cbind(diamstruct_forest1, Total)
write.table(diamstruct_forest, "E:\\Tellus\\ForMMI_forest_evaluation\\ref\\mean_diamstruct_forest.out", sep = "\t", row.names = FALSE, quote = FALSE)

