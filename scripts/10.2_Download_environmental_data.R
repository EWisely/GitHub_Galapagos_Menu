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

earlybiochem1of1<-stars::read_stars("../000_environmental_data/cmems_mod_glo_bgc_my_0.25_P1D-m_1718137929575.nc")
#chl, no3, nppv, o2, po4, si,
earlychl<-earlybiochem1of1$chl
str(earlychl)



earlybiochem<-terra::rast("../000_environmental_data/cmems_mod_glo_bgc_my_0.25_P1D-m_1718137929575.nc")



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


for (i in 1:nrow(pts)){pts$envlayer[i]=which(pts$unixdate[i]==biochem5of6@ptr$time)[1]}

pts<-pts %>% dplyr::mutate(envlayer1 = replace_na(envlayer, 1))


plot(earlybiochem["nppv"])
plot(earlybiochem["si"])

earlybiochem@ptr$time
which(pts$unixdate[1]==earlybiochem@ptr$time)[1]


#chl, no3, nppv, o2, po4, si,
early_chl<-earlybiochem["chl"]
late_chl<-biochem6of6["chl"]
early_no3<-earlybiochem["no3"]
late_no3<-biochem5of6["no3"]
early_nppv<-earlybiochem["nppv"]
late_nppv<-biochem3of6["nppv"]
early_o2<-earlybiochem["o2"]
late_o2<-biochem3of6["o2"]
early_po4<-earlybiochem["po4"]
late_po4<-biochem5of6["po4"]
early_si<-earlybiochem["si"]
late_si<-biochem5of6["si"]

for (i in 1:nrow(pts)){pts$early_envlayer[i]=which(pts$unixdate[i]==earlybiochem@ptr$time)[1]}


chltest = terra::extract(early_chl, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

chltest

no3test = terra::extract(early_no3, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

no3test

nppvtest = terra::extract(early_nppv, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

nppvtest

o2test = terra::extract(early_o2, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

o2test

po4test = terra::extract(early_po4, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

po4test

sitest = terra::extract(early_si, 
                         pts%>%
                           select(DD_long, DD_lat,)%>%
                           rename(lon=DD_long, lat=DD_lat),
                         method="bilinear",
                         layer=pts$early_envlayer[1])

sitest

### add the later data for 2022-2023
which(pts$unixdate[58]==late_si@ptr$time)[1]
pts$envlayer[58][1]


silate = terra::extract(late_si, 
                        pts%>%
                          select(DD_long, DD_lat)%>%
                          rename(lon=DD_long, lat=DD_lat),
                        method="bilinear",
                        layer=pts$envlayer1[1])

pts$envlayer[38][1]

silate

o2late = terra::extract(late_o2, 
                        pts%>%
                          select(DD_long, DD_lat)%>%
                          rename(lon=DD_long, lat=DD_lat),
                        method="bilinear",
                        layer=pts$envlayer1[1])


nppvlate = terra::extract(late_nppv, 
                        pts%>%
                          select(DD_long, DD_lat)%>%
                          rename(lon=DD_long, lat=DD_lat),
                        method="bilinear",
                        layer=pts$envlayer1[1])


po4late = terra::extract(late_po4, 
                          pts%>%
                            select(DD_long, DD_lat)%>%
                            rename(lon=DD_long, lat=DD_lat),
                          method="bilinear",
                          layer=pts$envlayer1[1])

no3late = terra::extract(late_no3, 
                          pts%>%
                            select(DD_long, DD_lat)%>%
                            rename(lon=DD_long, lat=DD_lat),
                          method="bilinear",
                          layer=pts$envlayer1[1])

chllate = terra::extract(late_chl, 
                          pts%>%
                            select(DD_long, DD_lat)%>%
                            rename(lon=DD_long, lat=DD_lat),
                          method="bilinear",
                          layer=pts$envlayer1[1])





#add to pts
pts<-cbind(pts, sitest$value)
pts<-cbind(pts, po4test$value)
pts<-cbind(pts, o2test$value, nppvtest$value)
pts<-cbind(pts, no3test$value, chltest$value)

pts<-cbind(pts, silate$value, po4late$value, o2late$value, nppvlate$value, no3late$value,chllate$value)

mutate(database=
         if_else(is.na(Scientific_name)==TRUE,
                 NA,
                 database))
pts1<-pts%>%
  mutate(silicate=
           if_else(envlayer1==1,
                   `sitest$value`,
                   `silate$value`),
         oxygen=
           if_else(envlayer1==1,
                   `o2test$value`,
                   `o2late$value`),
         phosphate=
           if_else(envlayer1==1,
                   `po4test$value`,
                   `po4late$value`),
         primary_productivity=
           if_else(envlayer1==1,
                   `nppvtest$value`,
                   `nppvlate$value`),
         nitrate=
           if_else(envlayer1==1,
                   `no3test$value`,
                   `no3late$value`),
         chloride=
           if_else(envlayer1==1,
                   `chltest$value`,
                   `chllate$value`))

pts2<-pts1%>%
  select(DD_long,DD_lat,Date.collected,sst,depth,silicate,oxygen,phosphate,primary_productivity,nitrate,chloride)

write.csv(pts2, file = "../000_environmental_data/Menu_sampling_envvars.csv")


