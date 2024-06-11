#Eldridge Wisely
#following tutorial on https://jcoliver.github.io/learn-r/011-species-distribution-models.html
#species distribution modeling for scalloped hammerheads from eDNA data in Galapagos shark nursery bays.

library(terra)
library(geodata)
library(predicts)
library(tidyverse)

# Read in hammerhead observations
Hammerhead_samples.ps <- readRDS("../10_Phyloseq/Menu_combined2/Hammerhead_samples_Phyloseq.rds")

obs_data<-as.data.frame(cbind(sample_data(Hammerhead_samples.ps)))
obs_data<-obs_data%>%
  select(DD_lat, DD_long, Date.collected)
  

# Check the data to make sure it loaded correctly
summary(obs_data)

# Store boundaries in a single extent object
#geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))
#for just the Galapagos: POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))
geographic_extent <- ext(x = c(-93.339844, -87.1875, -3.162456, 2.547988))

# Download data with geodata's world function to use for our base map
world_map <- world(resolution = 3,
                   path = "data/")

# Crop the map to our area of interest
my_map <- crop(x = world_map, y = geographic_extent)


# Make an extent that is 25% larger
sample_extent <- geographic_extent * 1.25

#grab bathymetry data
bathy <- terra::rast("../../../Physalia_GIS_in_R_2023/course_materials/practicals/data/bathymetry_galapagos/GEBCO_30_Oct_2023_672896985ee1/gebco_2023_n2.3839_s-2.9182_w-93.4343_e-85.3989.tif")


plot(bathy, col = hcl.colors(100, "blues"))

# Crop bathymetry data to desired extent
bathy_crop <- crop(x = bathy, y = sample_extent)


#send to PDF
pdf(file = "testplot.pdf")


# Plot the bathy
plot(bathy_crop, col = hcl.colors(100, "blues"))


# Plot the base map
plot(my_map,
     axes = TRUE, 
     col = "olivedrab", add=TRUE)



# Add the points for individual observations
points(x = obs_data$DD_long, 
       y = obs_data$DD_lat, 
       col = "red", 
       pch = 20, 
       cex = 1)


dev.off()






# Plot the first of the bioclim variables to check on cropping
plot(bathy_crop[[1]])



