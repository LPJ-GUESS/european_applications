######### Daily GPP modis to yearly GPP extract #########

library(raster)
library(rgdal)
library(terra)
library(DGVMTools)
library(rasterVis)

# Define new quantity
GUESS <- defineQuantity(
  id = "Europe",
  name = "Biomass",
  units = "kgC/m^2",
  format = "GUESS",
  add.to = GUESS
)
# Load LPJ-GUESS grid cells
output_dir1<-("C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\ref\\")
setwd(output_dir1)
ref_files <- list.files(output_dir1, pattern = "\\.out") # List the files in the output directory
print(ref_files)

# Establish a source object
EU_REF<- defineSource(id = "EU_reference",
                      dir = output_dir1, 
                      format = GUESS,
                      name = "Obs")

# Get gridcell field from source object
Grid.full <- getField(source = EU_REF, 
                      quant = "Europe")

# Plot check
Total_density <- plotSpatial(Grid.full, 
                             layers = c("Total"),
                             text.multiplier = 1.1, 
                             map.overlay = "world")
# Data and metadata check
print(Total_density)
print(Grid.full)
print(Grid.full@source)

# Get Lon Lat and project them
Lon<-Grid.full@data$Lon
Lat<-Grid.full@data$Lat
Lon_Lat<- data.frame(Lon, Lat)
coordinates(Lon_Lat) <- c("Lon", "Lat")
proj4string(Lon_Lat) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Make extent object from gridcells
tellus_ext<-Grid.full@spatial.extent

# Read LC and crop to established extent
corine <- raster("C:\\Users\\Admin\\Documents\\tellus\\corine\\u2018_clc2018_v2020_20u1_raster100m\\DATA\\U2018_CLC2018_V2020_20u1.tif")
corine <- projectRaster(corine, crs = "+proj=longlat +datum=WGS84 +no_defs", method = "ngb")
corine_crop1 <- crop(corine, tellus_ext)
writeRaster(corine_crop1,filename=file.path("C:\\Users\\Admin\\Documents\\tellus\\check_gpp\\cropped_corine.tif"), format = "GTiff", overwrite = TRUE)
mask_raster <- raster(res=res(corine_crop1), xmn=xmin(corine_crop1), xmx=xmax(corine_crop1), ymn=ymin(corine_crop1), ymx=ymax(corine_crop1), crs=crs(corine_crop1))

# 
# Filter the attribute table to only include the desired land use codes
mask_codes <- c(23, 24, 25) # Define the desired land use IDs
# 
#mask_raster <- raster(res=res(corine), xmn=xmin(corine), xmx=xmax(corine), ymn=ymin(corine), ymx=ymax(corine), crs=crs(corine))
values(mask_raster) <- 0

# Set the values to 1 where the land cover code is in the list of desired codes
for (i in mask_codes) {
  mask_raster[corine_crop1 == i] <- 1
}
plot(mask_raster)
writeRaster(mask_raster,filename=file.path("C:\\Users\\Admin\\Documents\\tellus\\check_gpp\\mask.tif"), format = "GTiff", overwrite = TRUE)

# Loading forest weights and crop to set extent
n<-raster("C:\\Users\\Admin\\Documents\\tellus\\check_gpp\\vikter.tif")
n<-crop(n,tellus_ext,snap="near")
n[n<0.8] <- 0
n[n>=0.8]<-1
#n<-projectRaster(n, crs = "+proj=longlat +ellps=sphere +no_defs", method = "ngb")
n@extent@xmax <- gpp_crop@extent@xmax
n@extent@xmin <- gpp_crop@extent@xmin
n@extent@ymax <- gpp_crop@extent@ymax
n@extent@ymin <- gpp_crop@extent@ymin

# Load yearly GPP files
gpp_files <- list.files("C:\\Users\\Admin\\Documents\\tellus\\GPP_GF\\", pattern=".tif$", full.names=TRUE)

#create an empty data table to store the final output
output_dt <- data.table(Lon = numeric(), Lat = numeric(), Year = integer(), Total = numeric())
nomask_dt <- data.table(Lon = numeric(), Lat = numeric(), Year = integer(), Total = numeric())
#nplot <- data.table(Lon = numeric(), Lat = numeric(), Year = integer(), Total = numeric())

# Loop over the annual GPP files, cropping, applying wights and aggregate. Creating a weighted annual, un-weighted annual & nplot output.
for (file in gpp_files) {
  year <- substr(file, nchar(file) - 7, nchar(file)-4)
  gpp_raster <- raster(file) #read GPP file
  gpp_crop <- crop(gpp_raster, tellus_ext)
  #gpp_crop <- projectRaster(gpp_crop,crs = "+proj=longlat +datum=WGS84 +no_defs",method = "ngb")
  #plot(gpp_crop/1000,main=paste("GPP no mask (0.05)",year))
  
  gpp_masked<-gpp_crop*n
  gpp_masked[gpp_masked==0]<-NA
  gpp_masked<-gpp_masked/1000
  gpp_crop1 <- gpp_crop/1000
  gpp_crop1[gpp_crop1==0]<-NA
  #plot(gpp_masked,main=paste("GPP(0.05)", year))
  gpp_nomask <- aggregate(gpp_crop1, fact=10,fun=mean, na.rm=T) 
  gpp_agg <- aggregate(gpp_masked, fact=10,fun=mean, na.rm=T)
  plot(gpp_nomask,main=paste("No Mask GPP (0.5)", year))
  plot(gpp_agg,main=paste("Mask GPP (0.5)", year))
  
  extract_values <- extract(gpp_agg, Lon_Lat, na.rm = TRUE) 
  gpp_crop_values <- extract(gpp_nomask, Lon_Lat, na.rm=TRUE)
  #nplot_values <- extract(nplots, Lon_Lat, na.rm = TRUE)
  
  #convert the grid list coordinates to data table format and add the year and GPP values
  year_vector <- rep(year, length(extract_values))
  output_dt_year <- data.table(x = Lon_Lat$Lon, y = Lon_Lat$Lat, Year = year_vector, Total = extract_values)
  #nplot_year <- data.table(x = Lon_Lat$Lon, y = Lon_Lat$Lat, Year = year_vector, Total = nplot_values)
  nomask_dt_years <-data.table(x = Lon_Lat$Lon, y = Lon_Lat$Lat, Year = year_vector, Total = gpp_crop_values )
  setnames(output_dt_year, c("Lon", "Lat", "Year", "Total"))
  setnames(nomask_dt_years, c("Lon", "Lat", "Year", "Total"))
  #setnames(nplot_year, c("Lon", "Lat", "Year", "Total"))
  nomask_dt <- rbind(nomask_dt,nomask_dt_years)
  #nplot <- rbind(nplot,nplot_year)
  output_dt <- rbind(output_dt, output_dt_year)
  
}

# Save outputs
write.table(output_dt, "C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\ref\\agpp_Europe_Torbern.out", sep = "\t",row.names = FALSE, quote = FALSE)
write.table(nplot, "E:\\Tellus\\ForMMI_forest_evaluation\\ref\\agpp_nplot.out", sep = "\t",row.names = FALSE, quote = FALSE)
write.table(nomask_dt, "C:\\Users\\Admin\\Documents\\tellus\\ForMMI_forest_evaluation\\ref\\agpp_Europe1.out", sep = "\t",row.names = FALSE, quote = FALSE)

gpp_masked1 <- raster(res=res(gpp_masked), xmn=xmin(gpp_masked), xmx=xmax(gpp_masked), ymn=ymin(gpp_masked), ymx=ymax(gpp_masked), crs=crs(gpp_masked))
values(gpp_masked1) <- 0
gpp_masked1[gpp_masked > 0] <- 1
gpp_masked1[gpp_masked1==0]<-NA
writeRaster(gpp_masked1,filename=file.path("E:\\Tellus\\check_gpp\\nplots.tif"), format = "GTiff", overwrite = TRUE)
nplots<-raster("E:\\Tellus\\check_gpp\\nplots.tif")
plot(nplots)
nplot_values <- extract(nplots, Lon_Lat, na.rm = TRUE)
