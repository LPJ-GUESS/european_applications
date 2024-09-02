# Script to extract lat lon, an approximated disturbance year and biome
# for TELLUS benchmarking runs with LPJ-GUESS.

# 02.12.2022
# Annemarie Eckes-Shephard

#Altered for Tellme europe regrowth simulations for ForestPaths/EuropeanApplications paper
# It now creates outputs files 
# 1) that follows this the LPJ-GUESS land use change structure
# 2) that is the gridlist file to force gegrowth simulations with.
# 14.08.2024- 26.08.2024
#Annemarie Eckes-Shephard

# No disturbance year records were found in the data,
# so to calculate the disturbance year
# I use the publication year - 5 years 
# (assuming time for data cleaning, obtaining and writing of publication) 
# as the reference year, from which I subtract the maximum age in the data.
# for Teobaldelli (Cannell) that was not possible (see respective notes below)

# Biome data for subsetting Teobaldelli dataset into temperate and boreal
# from : https://ecoregions.appspot.com/

# for the FIN data, disturbance years were created based on stand age and publication date of Repo et al 2019.
# gridcells to be disturbed were seleced randomply from the set of gridcells that cover
# mid and southern FIN in the standard lpjg gridlist for european application runs in 2024 Aug ( 0.5 * res)


####################################################################
library(readxl)
library(dplyr)
library(maps)
library(RColorBrewer)
path <- paste0(getwd(),"/data/")

#algorithm to turn site tropical lat longs into cru -gridcell lat lons
# based on translate:gridlist_to_coord in cruinput.cpp
#needs lat lon as "Lat" and "Lon"
#NB this function overwrites the values with updated "cru-compatible" values, to double-check that 
# some tropical gridcells don't end up in the same gridcell. If this is the case, I remove one of these data
# points.
create_crucompatible_gridlist <- function(df){
  soilinput.STEP = 0.5
  df$Lon = floor(df$Lon * (1/ soilinput.STEP)) * soilinput.STEP + soilinput.STEP / 2.0 -0.25
  df$Lat = floor(df$Lat * (1 / soilinput.STEP)) * soilinput.STEP + soilinput.STEP / 2.0 -0.25
  return(df)
}

####################################################################
#Temperate and boreal regrowth, Teobaldelli et al :
# biomass from here, as done in Pugh et al 2019, preprocessed by Vanessa Haverd
biom_data <- read.csv(paste0(path,"/Teobaldelli/Teobaldelli_biomassdata3.csv"))
# age from here, as done in Pugh et al 2019,, preprocessed by Vanessa Haverd
plot_info <- read.csv(paste0(path,"/Teobaldelli/Teobaldelli_plotinfo4.csv"))
# xlsx files - for merge to get species and country, used to select regions and relevant species for regrowth comparison:
Teo <- read_excel(paste0(path,"/Teobaldelli/Teobaldelli.xlsx"), sheet = "Plot Info")
# remove records which may only have root or leaf biomass records
# note: Removing these records is also important to get the woody biomass for the regrowth curve right
# remove "compartment 5" to only take into account aboveground woody biomass
biom_data <- biom_data[which(biom_data$Compartment..1.leaf.2.stem.branches.3.branches.4.stem.5.roots.6.b_bark.7.s_bark.8.other.!=5),]
#also rm leaf biomass, values are insignificant, but better for adequate comparison against model output:
biom_data <- biom_data[which(biom_data$Compartment..1.leaf.2.stem.branches.3.branches.4.stem.5.roots.6.b_bark.7.s_bark.8.other.!=1
& biom_data$Compartment..1.leaf.2.stem.branches.3.branches.4.stem.5.roots.6.b_bark.7.s_bark.8.other.!=5),]
#rm observations for which there is missing data for our relevant variables, biomass,Age:
plot_info <- plot_info[which(!is.na(plot_info$А..years.)),]
mer <- merge(biom_data,plot_info, by="ID")
# full merge
Teo_full_merge <- merge(mer,Teo, by="ID")
# new column to determine max age
Teo_full_merge_tmp <- Teo_full_merge %>% group_by(lat,lon) %>%  mutate(max_age = max(А..years.))#

#remove missing values
Teo_full_merge_tmp = Teo_full_merge_tmp[which(Teo_full_merge_tmp$lon > -6000), ]

#Shapefile from ecoregions
library(sf)
# from : https://ecoregions.appspot.com/
ecoregions <- st_read("/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs/Ecoregions2017/Ecoregions2017.shp")


# subset for biomes I want:
eco_Temperate <- ecoregions[which(ecoregions$BIOME_NUM == 4), ] 
                                 # |ecoregions$BIOME_NUM == 5), ] #  "Temperate Broadleaf & Mixed Forests" &  "Temperate Conifer Forests"

#transform coordinate system of TEo-dataset to the same as the ecoregions them to the system you want
sf_use_s2(FALSE)
Teo_full_merge_tmp_coords <- Teo_full_merge_tmp[!is.na(Teo_full_merge_tmp$lon & Teo_full_merge_tmp$lat),]
coords_Teo <- st_as_sf(Teo_full_merge_tmp_coords, coords = c("lon", "lat"), crs = 4326) 

# filter for places within temperate zone for conifers and broadleaves:
tmp_temperate <- coords_Teo %>% st_intersection(eco_Temperate)

Teo_full_merge_tmp_coords$Biome <- NA
Teo_full_merge_tmp_coords[(Teo_full_merge_tmp_coords$ID %in% tmp_temperate$ID),]$Biome <- "Temperate"

#remove datapoints with non-forested/non-Temperate biomes:
Teo_full_merge_tmp_coords <- Teo_full_merge_tmp_coords[which(!is.na(Teo_full_merge_tmp_coords$Biome)),]

#quick sanity checking where points cluster:
plot(Teo_full_merge_tmp_coords$lon,Teo_full_merge_tmp_coords$lat)
points(Teo_full_merge_tmp_coords[which(Teo_full_merge_tmp_coords$Biome=="Temperate"),]$lon,Teo_full_merge_tmp_coords[which(Teo_full_merge_tmp_coords$Biome=="Temperate"),]$lat, col="green")

#bring together (legacy from Tellus-Preparation code,where I have boreal data,too but done here to keep as similar as possible)
Teo_full_merge_tmp <- Teo_full_merge_tmp_coords

#################################

#Usoltsev 2001
#"The first summary materials on the forest biomass for the former USSR
#included 27 empirical defintions (Rodin Bazilevich 1965). After 1970 the fullest for that period for
#the above-mentioned territory were published in the amount of 185 definitions,
#then for the forests .. Utrkin 1970,.. Pozdnyakov 1975..
#To date the fullest database for biomass fractions ..."
#In the introduction it reads that the Usoltev dataset seems to be a new/updated database 
#for biomass fractions, which was then published after all the others.
#published in 2001, assuming 5 year assemblage period, assuming that latest data is from 1994
# and counting all sites backwards from that..
#initialise reference year column:
Teo_full_merge_tmp$refyear <- NA
Teo_full_merge_tmp[Teo_full_merge_tmp$Reference.nr==1,]$refyear <- 1994
#################################
#Cannell et al 1982 used main sources ( see introduction "scope" in
#https://daac.ornl.gov/daacdata/global_vegetation/biomass_allocation/comp/)
#Rodin an Bzxilevich 1976, Duvigneaud 1971, Young 1971,Ellenberg 1971
# Art and MArks 1971, Golley and Golley 1972,Young 1973,Golley and MEdina 1975
#Lieth and Wittaker 1975,Young 1975,Shidei and Kira 1977,Lamotte and Bourliere 1978
#Lieth 1978,Parte 1980, Reichle 1981
# I will be taking the mean of these dates, and -as above-
# subtract 5 years,- max stand age to get the disturbance year
refyear <-  floor(mean(c(1976,1971,1971,1971,1971,1972,1973,1975,1975,1975,1977,1978,1978,1980,1981)))
Teo_full_merge_tmp[Teo_full_merge_tmp$Reference.nr==2,]$refyear <-refyear
Teo_final <- Teo_full_merge_tmp %>% group_by(lat,lon) %>% mutate(dist_year = refyear - max_age)
#check final data:
unique(Teo_final[c("refyear","max_age","lat","lon","dist_year")])

# some data is recorded both in Usoltev and Cannell.Based on the above assumptions for designating disturbance timings,
# this leads to duplicate stands, i.e. to the same stand having two different disturbance timings assigned. 
# will be removed below.


####################################################################

# Create preliminary Gridlist df (contains duplicated rows) due to the nature of the original datastet, 
# as we only need the dist year entry as row for our gridlist, remove the duplicates:
gridlist <- unique(data.frame(lon= c(Teo_final$lon),
                              lat = c(Teo_final$lat), 
                              dist_year = floor(c(Teo_final$dist_year)),
                              biome = (Teo_final$Biome) )
)


####to merge with european coordinates: 
gridlist$lon <- gridlist$lon+0.25
gridlist$lat <- gridlist$lat+0.25
gridlist$in_benchmark_gridlist <- 1

#only temperate europe:
gridlist_temperate<- gridlist  

#remove swedish lat lons (conversation with Tom Pugh) - productivity gradient strongly visible in TreeMort data.
gridlist_temperate <- gridlist_temperate[which(gridlist_temperate$lat<54),]
#points(gridlist_temperate$lon,gridlist_temperate$lat,col="green",cex=0.2,pch=16)


#Now read in european gridlist: 
gridlist_europe <- read.table(paste0(path,"gridlist.txt"))
names(gridlist_europe) <- c("lon","lat")


whats_left <- merge(gridlist_europe,gridlist_temperate,by=c("lon","lat"))#,all.x=TRUE)



map(database = "world", fill = T, col = 'lightgray', xlim=c(-20,40.5),ylim=c(30,90))
points(gridlist_europe$lon,gridlist_europe$lat,cex=0.06, col="blue")
points(whats_left$lon,whats_left$lat,cex=0.7,col="magenta",pch=16)


#points(gridlist_temperate$lon,gridlist_temperate$lat,col="red",pch=16) 

legend(legend=c("gridlist centers","Temperate regrowth","Boreal regrowth"), fill=c("blue","magenta","darkgreen"),"bottomright",cex=0.1)

whats_left$biome <- NULL
whats_left$in_benchmark_gridlist <- NULL
gridlist_temperate <- whats_left
gridlist_temperate$Biome <- "Temperate"
rm(whats_left)



#############Finland simulations:
gridlist_europe <- read.table(paste0(path,"gridlist.txt"))
names(gridlist_europe) <- c("lon","lat")

#subset for fin data based on lat lon:
gridlist_fin <- gridlist_europe[which(gridlist_europe$lon>=21.25 & (gridlist_europe$lat>60.0 & gridlist_europe$lat<64.0)),]
#points(gridlist_fin$lon,gridlist_fin$lat,cex=0.06, col="green")

#add more data further north:
fin2 <- gridlist_europe[which(gridlist_europe$lon>=23.25 & (gridlist_europe$lat>60.0 & gridlist_europe$lat<=65.0)),]
#points(fin2$lon,fin2$lat,cex=0.06, col="purple")
#based on Fig 1, Repo et al 2019

# merge the two together
gridlist_fin <- merge(gridlist_fin, fin2,all=TRUE)
gridlist_fin$dist_year <- NA 



#read in regrowth data to determine the disturbance years:
Regr<- read.csv(paste0(path,"/Finland/regrowth_unmanaged_forests_last_table_biomass-in-bins-v3.csv"))
#Regr2 <-Regr %>% mutate(n_scaled = n/ sum(Regr$n)*156)# sum(Regr$n)
Regr$age_num_middle <- c(10,30,50,70,90,110,130,150,170,190,210)

Regr$dist_year <- NA
gridlist_fin$age_class <- NA
ncells <- dim(gridlist_fin)[1]

idxvalues.unique <- 1:ncells
set.seed(123)
# Create a sample without replacement (i.e. take the ball out and don't put it back in)
i=1
while(length(idxvalues.unique)>14){
  sample1 <- sample(x = idxvalues.unique,size = 14,replace = FALSE)
  gridlist_fin[sample1,]$dist_year <- 2019-5-Regr$age_num_middle[i]-floor(runif(14,0,1)*20-10)# randomise age a bit, within 20 years of bin size.
  gridlist_fin[sample1,]$age_class <- i 
  idxvalues.unique <- idxvalues.unique[!(idxvalues.unique %in% sample1)]
  i=i+1
}
# here, we make sure all age classes have at least 14 instances of regrowing stands occurring across finland.
# We are left with 5 age classes that were not allocated a gridcell using this method.
# That is ok, because in a test run I see that two corodinate pairs fail (at the coast), so I can reallocate these two age classes to
# two other coordinates that were originally unallocated.
# the other 3 are allocated to age class 11, 10 and 9 to yield the most datapoints agains observed age class bins in the comparison.


############
#Using Nitrogen deposition for (21.25,60.75)
#Grid cell not found in /data/benchmark_data/2015_12_14/ndep/GlobalNitrogenDeposition.bin

#Error: could not find/read data for (21.25,60.75), skipping it!
#  Using Nitrogen deposition for (21.25,62.25)
#    Grid cell not found in /data/benchmark_data/2015_12_14/ndep/GlobalNitrogenDeposition.bin

#Error: could not find/read data for (21.25,62.25), skipping it!
#  Using Nitrogen deposition for (21.25,62.75)
    

# coordinates not allocated any ages, at the coast
#points(24.25,64.25,col="red", cex=3)
#points(24.75,64.25,col="red", cex=3)
#points(24.75,64.75,col="red", cex=3)


# the first two will fail the LPJG-simulation because the gridlist coordinates are too close to the coast.
#Therefore, I reassign them to some leftover gridcells further inland, and north:
#ac_gc1 = gridlist_fin[which(gridlist_fin$lon==21.25 & gridlist_fin$lat ==60.75 ),]$age_class
#ac_gc2 = gridlist_fin[which(gridlist_fin$lon==21.25 & gridlist_fin$lat ==62.25 ),]$age_class
#gridlist_fin[idxvalues.unique[1],]$age_class <- ac_gc1
#gridlist_fin[idxvalues.unique[1],]$dist_year <- 2019-5-Regr$age_num_middle[ac_gc1]
#gridlist_fin[idxvalues.unique[2],]$age_class <- ac_gc2
#gridlist_fin[idxvalues.unique[2],]$dist_year <- 2019-5-Regr$age_num_middle[ac_gc2]

#handle the leftovers 
gridlist_fin[idxvalues.unique[1],]$age_class <- 7
gridlist_fin[idxvalues.unique[1],]$dist_year <- 2019-5-Regr$age_num_middle[7]
gridlist_fin[idxvalues.unique[2],]$age_class <- 8
gridlist_fin[idxvalues.unique[2],]$dist_year <- 2019-5-Regr$age_num_middle[8]
gridlist_fin[idxvalues.unique[3],]$age_class <- 11
gridlist_fin[idxvalues.unique[3],]$dist_year <- 2019-5-Regr$age_num_middle[11]
gridlist_fin[idxvalues.unique[4],]$age_class <- 10
gridlist_fin[idxvalues.unique[4],]$dist_year <- 2019-5-Regr$age_num_middle[10]
gridlist_fin[idxvalues.unique[5],]$age_class <- 9
gridlist_fin[idxvalues.unique[5],]$dist_year <- 2019-5-Regr$age_num_middle[9]

gridlist_fin$Biome <- "Boreal"

# Generate a palette of distinct colors (e.g., "Set3" has a good number of distinct colors)
colors <- brewer.pal(n = 11, name = "Set3")

# Assign colors to each age class
gridlist_fin$color <- colors[gridlist_fin$age_class]

#visual inspection of randomness in sampling of age classes in space.
map(database = "world", fill = T, col = 'lightgray', xlim=c(10.0,40.5),ylim=c(60.0,70.0))
points(gridlist_europe$lon,gridlist_europe$lat,cex=0.06, col="blue")
points(gridlist_fin$lon,gridlist_fin$lat,cex=0.7,col=gridlist_fin$color,pch=16)

legend(legend=c(as.character(sort(unique(gridlist_fin$age_class)))),col=c(colors[1:11]),pch=c(rep(16,11)),"topright")

#removing column that was only used for testing:
gridlist_fin$age_class <- NULL
gridlist_fin$color <- NULL
gridlist_out <- gridlist_fin

gridlist_out <- rbind(gridlist_fin,gridlist_temperate)


#last visual sanity-check before writing out:
png(filename = "data/processed/Appendix_regrowth.png",width=95,height=190,units = "mm",res=1260)

map(database = "world", fill = T, col = 'lightgray', xlim=c(-20,40.5),ylim=c(30,90))
points(gridlist_europe$lon,gridlist_europe$lat,cex=0.06, col="blue")
points(gridlist_fin$lon,gridlist_fin$lat,cex=0.4,col="darkgreen", pch=16)
points(gridlist_temperate$lon,gridlist_temperate$lat,cex=0.4,col="magenta",pch=16)
legend(legend=c("gridlist centers","Temperate regrowth","Boreal regrowth"), fill=c("blue","magenta","darkgreen"),"bottomright",cex=0.5)
dev.off()


#create output file that follows this the LPJ-GUESS land use change structure:
#Year	NATURAL	FOREST
#1700	1	0
#1701	0	1
disturbance_year <- gridlist_out
disturbance_year$FOREST  <- 1
#disturbance_year$FOREST  <- 0
disturbance_year$NATURAL <- 0
disturbance_year$BARREN <- 0

regrowth_year         <- gridlist_out
regrowth_year$FOREST  <- 0
regrowth_year$NATURAL <- 1
regrowth_year$BARREN <- 0

#adressing error "Make sure all coordinates have the same data years !";
# in lpjg:
# Create a complete data frame for all lon, lat, and year combinations
complete_cases <- regrowth_year %>%
  # Group by lat and lon to focus on existing combinations
  group_by(lon, lat) %>%
  # Define all possible combinations of lon, lat, and year
  complete( dist_year = seq(from = min(regrowth_year$dist_year), to = 2023)) %>%
  # Fill in FOREST and NATURAL with the last observed value for each lon and lat
  
  fill(FOREST, NATURAL, .direction = "downup") %>%
  ungroup() %>%
  # Replace any NA values in FOREST with 0 and in NATURAL with 1
  mutate(
    FOREST = ifelse(is.na(FOREST), 0, FOREST),
    NATURAL = ifelse(is.na(NATURAL), 1, NATURAL),
    BARREN = ifelse(is.na(BARREN), 0, BARREN)
  )


luc_tmp <- rbind(complete_cases,disturbance_year)
luc_tmp$Biome <- NULL

#find duplicates, remove:
filtered_data <- luc_tmp %>%
  group_by(lon, lat, dist_year) %>%
  filter(!(FOREST == 0 & NATURAL == 1 & 
             any(FOREST == 1 & NATURAL == 0))) %>%
  ungroup()

#sanity check
sum(duplicated(filtered_data[c("lat", "lon", "dist_year")]))
# is good

luc_out <- filtered_data %>% arrange(desc(lat),desc(lon),dist_year,desc(FOREST))
names(luc_out)[3]<- "year"


#sanity check
map(database = "world", fill = T, col = 'lightgray', xlim=c(0.0,40.5),ylim=c(50,90))
uq <-unique(luc_out[c("lat","lon")])
points(uq$lon,uq$lat)
#shoudl only cover finland and other european contries somewhat- good


#sanity check
plot(type="l",luc_out[which(luc_out$lat == as.numeric(uq[1,1]) & luc_out$lon == as.numeric(uq[1,2])), ]$year,luc_out[which(luc_out$lat == as.numeric(uq[1,1]) & luc_out$lon == as.numeric(uq[1,2])), ]$FOREST)
for(i in 2:dim(luc_out)[2]){
  #plot(luc_out[which(luc_out$lat == 61.75 & luc_out$lon == 29.75), ]$year,luc_out[which(luc_out$lat == 61.75 & luc_out$lon == 29.75), ]$FOREST)
  lines(luc_out[which(luc_out$lat == as.numeric(uq[i,1]) & luc_out$lon == as.numeric(uq[i,2])), ]$year,luc_out[which(luc_out$lat == as.numeric(uq[i,1]) & luc_out$lon == as.numeric(uq[i,2])), ]$FOREST, col=i)
}
    
#some stats for the manuscript text:
dim(gridlist_out[which(gridlist_out$Biome=="Temperate"),])
min(gridlist_out[which(gridlist_out$Biome=="Temperate"),]$dist_year)
max(gridlist_out[which(gridlist_out$Biome=="Temperate"),]$dist_year)

dim(gridlist_out[which(gridlist_out$Biome=="Boreal"),])
min(gridlist_out[which(gridlist_out$Biome=="Boreal"),]$dist_year)
max(gridlist_out[which(gridlist_out$Biome=="Boreal"),]$dist_year)

#Create input files for running LPJ-GUESS:
write.table(luc_out,file="data/processed/luc_disturbance26082024.txt",
            col.names = TRUE, row.names =FALSE)

write.table(gridlist_out,file="data/processed_for_LPJG/gridlist_update29082024.txt",
            col.names = FALSE, row.names =FALSE)

##################






