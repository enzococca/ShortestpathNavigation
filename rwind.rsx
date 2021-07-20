##load_vector_using_rgdal
##load_raster_using_rgdal
##rwind=group
##Navigazione_Antica=name
##Anno=number 2020
##Mese=number 1
##Giorno=number 1
##Ora=number 0
##lon1=number -8
##lon2=number 38
##lat1=number 12
##lat2=number 50
#Hfactor=file
##Points_A=vector point
##Points_B=vector point
##Mask=vector polygon

#Labels_A=string A
#Labels_B=string B
#Output_speed=output raster
#Output_dist=output raster
##Output_A_to_B=output vector
#Output_B_to_A=output vector
##Output_speed_direction_vector=output vector
#Output_direction_vector=output vector
#showplots
library(sp)
library(sf)
library(maptools)
library(rWind)
library(raster)
library(gdistance)
library(rworldmap)

w <- wind.dl(Anno, Mese, Giorno, Ora, lon1, lon2, lat1, lat2)
w
p1 <- wind2raster(w)
#p1.gcs<-projectRaster(p1, crs="+init=EPSG:3857")

#projectRaster(p1, crs="+init=EPSG:3034")
#b <- as(extent(Mask), 'SpatialPolygons')
p <-(wind_layer <- do.call(stack,lapply(1:nrow(Mask), function(x) 
raster::mask(crop(p1, extent(Mask)), Mask[x,]) )) )
p
#t<-read.delim(Hfactor, header=FALSE, sep = "\t", dec = ".")

cd<-flow.dispersion(p,type="actvie",output="transitionLaye")
Conductance<- geoCorrection(cd, type="c")
AtoB<- shortestPath(Conductance,Points_A, Points_B, Cost, output="SpatialLines")
BtoA<- shortestPath(Conductance,Points_B, Points_A, output="SpatialLines")

raster.sp <- as(p$speed, "SpatialPixelsDataFrame") 

#Output_speed=raster.sp

raster2.sp <- as(p$direction, "SpatialPixelsDataFrame") 
#crs(raster2.sp)<-crs(p)
#Output_dist=raster2.sp
x <- SpatialLinesDataFrame(AtoB, data.frame(row.names=row.names(AtoB),
                        ID=1:length(AtoB)))
  class(AtoB)
y <- SpatialLinesDataFrame(BtoA, data.frame(row.names=row.names(BtoA),
                        ID=1:length(BtoA)))
  class(BtoA)

a1.sp<-as(x, "SpatialLinesDataFrame")

Output_A_to_B=a1.sp

b1.sp<-as(y, "SpatialLinesDataFrame")
#Output_B_to_A=b1.sp
spe<-rasterToPoints(p, spatial=TRUE)
#crs(spe)<-crs(p)
#spe<-SpatialPoints(spe)
#class(spe)
#dir<-rasterToPoints(p$direction, spatial=TRUE)
#crs(dir)<-crs(p)
#dir
#spp<-over(spe,raster2.sp)
Output_speed_direction_vector=spe
#Output_direction_vector=dir