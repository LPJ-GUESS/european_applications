# Script to map European species (tree and shrub) to Plant Functional Type

# 05.04.2024
# Haoming Zhong


#############################################################################
library(openxlsx)
library(dplyr)
library(stringr)
library(tidyverse)


# Read NFI data from each EU country
# Select species column and merge them together 
# Create a table called ** specieslist.xlsx **
# Based on species chorological maps and experts opinion
# Indentified species as Boreal and Temperate species


###### Process species list
Species <- read.xlsx('FinalSpecieslist.xlsx')
# Filter total fit species - the species already has a PFT, fx, Abies alba : Abi_alb.
Totalfit <- read.xlsx('totalfitlist.xlsx') 
Species<- anti_join(Species,Totalfit, by = "Species")

# Read reference shade tolerance data
ST <- read.csv('ShadetoleranceNiinemets.csv',sep=';')
# Read reference wood density data
WD <- read.delim('WD.txt', sep = '\t')
#Read reference leaf data
LE <- read.csv('TRY_categorical_traits_with_tolerances_Stefan.csv',sep=',')

tree <- Species %>%
  left_join(LE, by = c('Species' = 'AccSpeciesName')) %>%     
  left_join(ST, by = 'Species') %>% 
  left_join(WD, by = c('Species' = 'species'))
  
# select useful columns
tree <- tree[, c("Species", "Type",'Genus','Family','LeafType','LeafPhenology','Shade.tolerance.y','WD')]

write.table(tree,'treewithmissingdata.xlsx')

# There are some missing traits value, based on other reference such as articles, AI, or average Genus value
# We got a new comprehensive data set. 

treenew <- read.xlsx('Specieswithallinfo.xlsx')

####### Seperate Boreal and Temperate Species

borealsp <- treenew[treenew$Type=='Boreal',]
temperatesp <- treenew[treenew$Type=='Temperate',]


###### Boreal Species ######################

# Mapping Boreal species to Boreal PFTs
BorealPFT <- read.xlsx('BorealPFT.xlsx')

# Add value range to Wood density and Shade Tolerance (±10%)
# Wood density
BorealPFT$WD <- as.numeric(BorealPFT$WD)
BorealPFT <- BorealPFT %>%
  mutate(WD_lower = WD * 0.9,  
         WD_upper = WD * 1.1)

#Shade tolerance
BorealPFT <- BorealPFT %>%
  mutate(
    ShadeTolerance_values = str_split(ShadeTolerance, "-"),  
    ShadeTolerance_lower = as.numeric(sapply(ShadeTolerance_values, function(x) as.numeric(x[1]))),  
    ShadeTolerance_upper = as.numeric(sapply(ShadeTolerance_values, function(x) as.numeric(x[2])))   
  )
BorealPFT <- BorealPFT  %>%
  select(-ShadeTolerance_values)


####### Step1: match Leaf 
#########################
Boreal1 <- borealsp %>%
  inner_join(BorealPFT, by = c("LeafType", "LeafPhenology"), relationship = "many-to-many")

# Find only one species - one PFT
BorealResult1 <- Boreal1  %>%
  group_by(Species) %>%
  filter(n() == 1)  

# Delete the species already had one PFT
Boreal2 <- Boreal1 %>%
  anti_join(BorealResult1, by = "Species")  

# Check any missing species has no leaf match 
species_original1 <- unique(borealsp$Species)
species_aftermatch1 <- unique(Boreal1$Species)
Bspecies_missing1 <- as.data.frame(setdiff(species_original1, species_aftermatch1)). # 2 species with no leaf match

####### Step2: match Shade Tolerance
####################################
# filter rows by st, if within PFT's range then keep it, otherwise delete row
Boreal3 <- Boreal2 %>%
  filter(ShadeTolerance.x >= ShadeTolerance_lower & ShadeTolerance.x <= ShadeTolerance_upper)

# Find only one species - one PFT
BorealResult2 <- Boreal3  %>%
  group_by(Species) %>%
  filter(n() == 1)  # 只保留每个物种只有一行的数据

# Delete the species already had one PFT
Boreal4 <- Boreal3 %>%
  anti_join(BorealResult2, by = "Species")  

# Check any missing species has no ST match
species_original2 <- unique(Boreal2$Species)
species_aftermatch2 <- unique(Boreal3$Species)
Bspecies_missing2 <- as.data.frame(setdiff(species_original2, species_aftermatch2)) # 4 species with no ST match

####### Step3: match Wood Density
####################################
Boreal5 <- Boreal4 %>%
  filter(WoodDensity >= WD_lower & WoodDensity <= WD_upper)

BorealResult3 <- Boreal5  %>%
  group_by(Species) %>%
  filter(n() == 1)  

# Check any missing species has no WD match
species_original3 <- unique(Boreal4$Species)
species_aftermatch3 <- unique(Boreal5$Species)
Bspecies_missing3 <- as.data.frame(setdiff(species_original3, species_aftermatch3)) # no species with no WD match

# Write Boreal PFT result 1
borealPFTresult <- bind_rows(BorealResult1,BorealResult2,BorealResult3) # after 3-step matching, 28 species found the only one matched PFT

####### Step4: process missing species 
####################################

# Bspecies_missing1: 2 species with no leaf match, decision made by experts 

# Bspecies_missing2: 4 species with no ST match, check if they have WD match and Genus name match
BBspecies_missing2 <- merge(Boreal1,Bspecies_missing2,by.x='Species',by.y='setdiff(species_original2, species_aftermatch2)',all.x = F)


Boreal6 <- BBspecies_missing2 %>%
  filter(WoodDensity >= WD_lower & WoodDensity  <= WD_upper)

BorealResult4 <- Boreal6   %>%
  group_by(Species) %>%
  filter(n() == 1)  # even though 4 species has no ST match, but they have WD match. 

# Write Boreal PFT result 2
borealPFTresult <- bind_rows(borealPFTresult,BorealResult4) # after Step 4, 4 species found only one matched PFT. 
write.csv(borealPFTresult,'BorealspeciesPFT.csv')

# Write 2 unsure Boreal species, decision made by experts
write.csv(Bspecies_missing1,'borealunsure.csv')





###### Temperate Species ######################
# Mapping temperate species to temperate PFT

TemperatePFT <- read.xlsx('TemperatePFT.xlsx')

# Add value range to Wood density and Shade Tolerance (±10%)
# Wood density
TemperatePFT$WD <- as.numeric(TemperatePFT$WD)
TemperatePFT <- TemperatePFT %>%
  mutate(WD_lower = WD * 0.9,  
         WD_upper = WD * 1.1)

# Shade tolerance
TemperatePFT <- TemperatePFT %>%
  mutate(
    ShadeTolerance_values = str_split(ShadeTolerance, "-"),  
    ShadeTolerance_lower = as.numeric(sapply(ShadeTolerance_values, function(x) as.numeric(x[1]))), 
    ShadeTolerance_upper = as.numeric(sapply(ShadeTolerance_values, function(x) as.numeric(x[2])))   
  )
TemperatePFT <- TemperatePFT  %>%
  select(-ShadeTolerance_values)


####### Step1: match Leaf 
#########################
Temperate1 <- temperatesp %>%
  inner_join(TemperatePFT, by = c("LeafType", "LeafPhenology"), relationship = "many-to-many")

# Find only one species - one PFT
TemperateResult1 <- Temperate1  %>%
  group_by(Species) %>%
  filter(n() == 1)    # no species has only one leaf match

# Check any missing species has no leaf match 
species_original1 <- unique(temperatesp$Species)
species_aftermatch1 <- unique(Temperate1$Species)
Tspecies_missing1 <- as.data.frame(setdiff(species_original1, species_aftermatch1)) # 2 species with no leaf match



####### Step2: match Shade Tolerance
####################################

Temperate2 <- Temperate1 %>%
  filter(ShadeTolerance.x >= ShadeTolerance_lower & ShadeTolerance.x <= ShadeTolerance_upper)

TemperateResult2 <- Temperate2  %>%
  group_by(Species) %>%
  filter(n() == 1)   # 22 species 

Temperate3 <- Temperate2 %>%
  anti_join(TemperateResult2, by = "Species")  

# Check any missing species has no ST match 
species_original2 <- unique(Temperate1$Species)
species_aftermatch2 <- unique(Temperate2$Species)
Tspecies_missing2 <- as.data.frame(setdiff(species_original2, species_aftermatch2)) # 29 species with no ST match

####### Step3: match Wood Density
####################################
unique(Temperate4$Species)
Temperate4 <- Temperate3 %>%
  filter(WoodDensity >= WD_lower & WoodDensity <= WD_upper)

TemperateResult3 <- Temperate4  %>%
  group_by(Species) %>%
  filter(n() == 1)  # 38 species

Temperate5 <- Temperate4 %>%
  anti_join(TemperateResult3, by = "Species")  

# Check any missing species has WD leaf match 
species_original3 <- unique(Temperate3$Species)
species_aftermatch3 <- unique(Temperate4$Species)
Tspecies_missing3 <- as.data.frame(setdiff(species_original3, species_aftermatch3)) # 33 species with no WD match


####### Step4: match Genus name
####################################
TemperateResult4 <-  Temperate5 %>%
  group_by(Species) %>%
  filter(Genus.name == genus.cor) # 4 species 

# 33 species with multiple Leaf, WD and ST PFT matches 
Temperate6 <- Temperate5 %>%
  anti_join(TemperateResult4, by = "Species")  


####### Step5: process missing species 
####################################

# Tspecies_missing3: 33 species with no WD match but with ST and Leaf match, process Genus name match now
#########################################################################################################
# get species values 
missing3 <- merge(Temperate1,Tspecies_missing3,by.x='Species',by.y='setdiff(species_original3, species_aftermatch3)',all.x = F)

# Genus name match
TemperateResult5 <-  missing3 %>%
  group_by(Species) %>%
  filter(Genus.name == genus.cor) # 4 species

# 29 species have multiple ST and Leaf match but with no WD and Genus name match 
Temperate7 <- missing3 %>%
  anti_join(TemperateResult5, by = "Species")  


# Tspecies_missing2: 29 species with no ST match but with Leaf match, process WD and Genus name match now
#########################################################################################################
# get species values
missing2 <- merge(Temperate1,Tspecies_missing2,by.x='Species',by.y='setdiff(species_original2, species_aftermatch2)',all.x = F)

# Wood density match 
missing2.1 <- missing2 %>%
  filter(WoodDensity >= WD_lower & WoodDensity <= WD_upper)

TemperateResult6<- missing2.1  %>%
  group_by(Species) %>%
  filter(n() == 1)  # 14 species

# Genus name match
Temperate8 <- missing2.1 %>%
  anti_join(TemperateResult6, by = "Species")  

TemperateResult7 <-  Temperate8 %>%
  group_by(Species) %>%
  filter(Genus.name == genus.cor) # 2 species

# 13 species have multiple leaf match but no ST WD and Genus name match
Temperate9 <- missing2.1 %>%
  anti_join(TemperateResult6, by = "Species")  

Temperate9 <- Temperate9 %>%
  anti_join(TemperateResult7, by = "Species")  # 1 species 

# 12 species has no ST match and no WD match, but only Leaf match, Check if they have Genus name match
species_original4 <- unique(missing2$Species)
species_aftermatch4 <- unique(missing2.1$Species)
species_missingmissing1 <- as.data.frame(setdiff(species_original4, species_aftermatch4))

missingmissing <- merge(Temperate1,species_missingmissing1,by.x='Species',by.y='setdiff(species_original4, species_aftermatch4)',all.x = F)

TemperateResult8 <-  missingmissing %>%
  group_by(Species) %>%
  filter(Genus.name == genus.cor) # 4 species 

Temperate10 <- missingmissing %>%
  anti_join(TemperateResult8, by = "Species")  # 8 species


# Write Temperate PFT result 
temperatePFTresult <- bind_rows(TemperateResult1, TemperateResult2, TemperateResult3, TemperateResult4,TemperateResult5,TemperateResult6,TemperateResult7,TemperateResult8)
write.csv(temperatePFTresult,'TemperatespeciesPFT.csv')

# Write 73 unsure Temperate species, decision made by experts
temperatePFTresultunsure <- bind_rows(Temperate6,Temperate7,Temperate9,Temperate10)

# some process to add two missing species with no Leaf match
species2 <- merge(temperatesp,Tspecies_missing1,by.x='Species',by.y='setdiff(species_original1, species_aftermatch1)',all.x = F)
missing_cols <- setdiff(colnames(temperatePFTresultunsure), colnames(species2))
species2[missing_cols] <- NA
species2 <- species2[, colnames(temperatePFTresultunsure)]
temperatePFTresultunsure <- rbind(temperatePFTresultunsure, species2)

write.csv(temperatePFTresultunsure,'temperateunsure.csv')


# The final species-PFT result list is included in the appendix of the paper.
# Table A2.1. LPJ-GUESS PFT characteristics and traits of the associated species