#Eldridge Wisely
#June 10, 2024
#Get bathymetry data from marmap and environmental variables from Copernicus and add them to the sample data of the phyloseq object.

#library(CopernicusMarine)
#library(marmap)
library(tidyverse)
library(phyloseq)
#library(stars)
library(terra)
library(geodata)
library(basemaps)
library(ggrepel)
library(rnaturalearth)
#library(tidyterra)
#library(xts)

biochem1of6 <- terra::rast("../000_environmental_data/cmems_mod_glo_bgc_anfc_0.25deg_static_1718047948686.nc")
#> deptho, deptho_lev, (bathymetry depth under geoid, model level of depth)
biochem2of6<- terra::rast("../000_environmental_data/cmems_mod_glo_bgc_anfc_0.25deg_static_1718047985502.nc")
#e1t, e2t, (cell dimension along x axis, then cell dimension along y axis), not sure I'll need these.
biochem3of6<- terra::rast("../000_environmental_data/cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m_1718047809206.nc")
#nppv, o2 (Total Primary Production of Phyto net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water, dissolved oxygen)
biochem4of6<-terra::rast("../000_environmental_data/cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m_1718047754996.nc")
#dissic, ph, talk, (dissolved inorganic carbon, pH, total alkalinity)
biochem5of6<-terra::rast("../000_environmental_data/cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m_1718047875675.nc")
#fe, no3, po4, si, (iron, nitrate, phosphate, silicate)
biochem6of6<-terra::rast("../000_environmental_data/cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m_1718047896152.nc")
#chl, phyc, (chlorophyll phytoplankton)

plot(biochem1of6["deptho"], col = hcl.colors(100), axes = TRUE)

earlybiochem1of1<-terra::rast("../000_environmental_data/cmems_mod_glo_bgc_my_0.25_P1D-m_1718137929575.nc")



#SSTmap<-stars::read_stars("../000_environmental_data/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1718054723180.nc")
#analysed_sst, analysis_error


SSTmap<-terra::rast("../000_environmental_data/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1718054723180.nc")

SSTmap
plot(SSTmap["analysed_sst"])

SSTmap@ptr$time

#https://www.epochconverter.com/

# timeInfo(SSTmap)
# 
# terra::values(SSTmap)
# 
# terra::varnames(SSTmap)

# # tidyterra
# glimpse(SSTmap, width = NULL, n = 10, max_extra_cols = 20)

#following answer from https://gis.stackexchange.com/questions/460747/data-extraction-for-time-series-transform-multiple-netcdf-files-to-spatraster-w

SSTmap<-SSTmap["analysed_sst"]
SSTmap
terra::time(SSTmap) |> head(10)



# # create random point for extraction
# #p <- sf::st_point(c(6.1, 50.5)) |> 
#   sf::st_sfc(crs = "epsg:4326") |> 
#   terra::vect()
# 
# p<- sf::st_point(c(	
#   -89.46920,-0.800391))|>
#   sf::st_sfc(crs = "epsg:4326") |>
#   terra::vect()
# 
# # extract values from raster stack based on feature input
# #values <- terra::extract(r, p, ID = FALSE) |> as.double()
# 
# values<- terra::extract(SSTmap, p, ID=FALSE)|> as.double()
# values
# 
# # get timestamps
# datetimes <- terra::time(SSTmap) |> as.POSIXct(origin = "1970-01-01", tz = "UTC")
# 
# # create xts object
# xts <- xts::xts(values, order.by = datetimes)
# 
# # inspect
# plot(xts, type = "h", col = "blue")
# 
# xts
# 
# ############

# library(ncdf4)
# obsdata = open.ncdf("obs1.nc")
# obsdatadates = as.Date(obsdata$dim$time$vals,origin = '1950-01-01')
# #Get the whole data first
# obsoutput = get.var.ncdf(obsdata, varid = 'tasmin')
# #Prepare your points of interest
# points_of_interest = data.frame(lat=seq(1,8,by=2),lon=c(1,5,3,6))
# #Subset your data accordingly
# data_at_point = apply(points_of_interest,1,function(x)obsoutput[x[1],x[2],])
# #Turn it into a dataframe
# data_at_point = as.data.frame(data_at_point)
# #Add the dates to the dataframe
# data_at_point$Date = obsdatadates



###########
Shark_Menu_taxa_merge_species_named_traits.hell.ps<-readRDS(file="../10_Phyloseq/Menu_combined2/Menu_combined2_hellinger_transformed_traits_added_taxa_named_Phyloseq.rds")


pts<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))
pts<-pts%>%
  select(DD_lat, DD_long, Date.collected)%>%
  dplyr::mutate(unixdate= as.numeric(as.POSIXct(Date.collected, format="%m/%d/%y",
                                                tz="GMT")))

pts

pts1<-pts%>%
  select(DD_long,DD_lat)%>%
  rename(lon=DD_long, lat=DD_lat)

pts1

v <- vect(cbind(pts1), crs="+proj=longlat")
radius <- terra::buffer(v, 10)



#for each observation:
#get the appropriate raster layer for each unixdate
which(pts$unixdate[59]==SSTmap@ptr$time)[1]

SSTmap[[766]]

for (i in 1:nrow(pts)){pts$sstlayer[i]=which(pts$unixdate[i]==SSTmap@ptr$time)[1]}



#extract environmental data at the lat long from the appropriate layer
ssttest = terra::extract(SSTmap, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$sstlayer)

ssttest
class(ssttest)

#add SST to pts

#convert Kelvin to Celcius
sstemps<-ssttest%>%
  dplyr::mutate(sst=value-273.15)

#put sst into pts
pts<-cbind(pts, sstemps$sst)



  

# #plan B average across layers
# SSTmap_mean<-terra::mean(SSTmap, na.rm = T)
# 
# plot(SSTmap_mean)
# 
# sstdata = terra::extract(SSTmap_mean, 
#                          pts%>%
#                            select(DD_long, DD_lat)%>%
#                            rename(lon=DD_long, lat=DD_lat))
# 
# 
# sstdata



# # convert POSIXct time to character, to please ggplot's facet_wrap()
# z1 = st_set_dimensions(z, 3, values = as.character(st_get_dimension_values(z, 3)))
# library(ggplot2)
# library(viridis)
## Loading required package: viridisLite


#Bathymetry data and map

#Load bathymetry raster for the whole Galapagos region downloaded from GEBCO: https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2019/grid_terms_of_use.html

bathy <- terra::rast("../../../Physalia_GIS_in_R_2023/course_materials/practicals/data/bathymetry_galapagos/GEBCO_30_Oct_2023_672896985ee1/gebco_2023_n2.3839_s-2.9182_w-93.4343_e-85.3989.tif")
plot(bathy, col = hcl.colors(100, "blues"))
plot(radius, add=TRUE)
bathy

#get depths under points

depths<-terra::extract(bathy, 
                                 pts%>%
                                   select(DD_long, DD_lat)%>%
                                   rename(lon=DD_long, lat=DD_lat),
                                 method="bilinear",
                                  touches=TRUE)
#add to pts
pts<-cbind(pts, depths$`gebco_2023_n2.3839_s-2.9182_w-93.4343_e-85.3989`)

pts<-pts%>%
  dplyr::rename(depth=`depths$\`gebco_2023_n2.3839_s-2.9182_w-93.4343_e-85.3989\``, sst=`sstemps$sst`)

#Add primary productivity
plot(biochem3of6["nppv"])
biochem3of6@ptr$time

for (i in 1:nrow(pts)){pts$envlayer[i]=which(pts$unixdate[i]==biochem3of6@ptr$time)[1]}

