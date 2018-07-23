#Uses data collated & summarised in PerfectEnemyofGood_mergeData.r

library(sf)
library(ggplot2)
#library(gridExtra)#for grid.arrange
#library(lme4) #for lmer

#options(stringsAsFactors=FALSE) # turn off automatic factor coersion

wd <- "D:/Box Sync/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)

###EDITME

currdate <- "170612_static_MZregion_2008_2012"

########################
#SET UP DATA
########################
#read data
#allCountiesDF <- read.csv("2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static.csv", header=TRUE)
allCountiesDF <- read.csv("2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static_plus_sgmzs.csv", header=TRUE)
#set categorical variables
allCountiesDF$Year <- as.factor(allCountiesDF$Year)
allCountiesDF$ADMIN_FIPS <- as.factor(allCountiesDF$ADMIN_FIPS)
allCountiesDF$Region <- as.factor(allCountiesDF$Region)
names(allCountiesDF)[which(names(allCountiesDF)=="Acres.of.Irrigated.Harvested.Cropland.as.Percent.of.All.Harvested.Cropland.Acreage...2012")] <- "PercentLandIrrigated"
#fill blanks
allCountiesDF$PercentLandIrrigated[is.na(allCountiesDF$PercentLandIrrigated)] <- 0 #all land
allCountiesDF$PercentCroplandthatisIrrigated[is.na(allCountiesDF$PercentCroplandthatisIrrigated)] <- 0 #just cropland
#subset out only sage grouse counties & conversion on land capability class 1 to6	
sgCountiesDF1to6 <- subset(allCountiesDF, SageGrouseCounty==1 & LCC=="LCC1to6")
#drop new years
sgCountiesDF1to6 <- subset(sgCountiesDF1to6, Year %in% 2008:2012)	
sgCountiesDF1to6 <- droplevels(sgCountiesDF1to6)	
#drop rows with NAs in rent & conversion prop
sgCountiesDF1to6  <- sgCountiesDF1to6 [complete.cases(sgCountiesDF1to6 [,c("NASS_LandRent_NonIrrigatedCropland","NASS_LandRent_RangelandNEG", "ConversionPropCropRangeV2logZIF001")]),]
sgCountiesDF1to6 <- droplevels(sgCountiesDF1to6)

###Model
mod10a <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + Popn*AREA_URBAN_ha) )

###Load shapefile of counties
counties <- read_sf("1 Inputs/Boundaries/SGCounties/SGCounties_overlapping_PACSandBreeding_NAD83alb.shp")
centroidlatlon <- st_centroid(counties) %>% st_coordinates()
counties$lat <- centroidlatlon[,2]
counties$lon <- centroidlatlon[,1]

###Merge with model residuals
names(mod10a$model)[4] <- "Year"
model_res <- cbind(sgCountiesDF1to6[,1:4], mod10a$residuals)

###Plot each year
latlonDF <- lapply(2009:2012, function(curryr){
  yr_res <- model_res[model_res$Year==curryr, ]
  yr_res_sf <- merge(counties, yr_res, by.x="fips", by.y="ADMIN_FIPS", all.x=TRUE)
  png(sprintf("3 Model output/spatial autocorrelation_%s.png", curryr))
  plot(yr_res_sf["mod10a.residuals"], main=sprintf("%s residuals", curryr))
  dev.off()
 return(yr_res_sf)[,c("ADMIN_FIPS", "Year", "mod10a.residuals", "lat", "lon")]
})

###Plot by latitude & longitude
latlonDF <- do.call(rbind, latlonDF)
png("3 Model output/spatial autocorrelation/Model_residuals_latlon.png")
par(mfrow=c(2,1))
with(latlonDF, plot(mod10a.residuals ~ lat))
with(latlonDF, plot(mod10a.residuals ~ lon))
dev.off()

###Model residuals to check
sink("3 Model output/spatial autocorrelation/Models_of_residuals_latlon.txt")
modlat <- lm(mod10a.residuals ~ lat + Year, data=latlonDF)
summary(modlat)
modlon <- lm(mod10a.residuals ~ lon + Year, data=latlonDF)
summary(modlon)
modlatlon <- lm(mod10a.residuals ~ lat + lon + Year, data=latlonDF)
summary(modlatlon)
modlatloni <- lm(mod10a.residuals ~ lat + lon + lat*lon + Year, data=latlonDF)
summary(modlatloni)
modlatlonib <- lm(mod10a.residuals ~ lat + lon + lat*lon , data=latlonDF)
summary(modlatlonib)
sink()