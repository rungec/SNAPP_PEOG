#This script summarises and combines dataset for the PERFECT ENEMY OF THE GOOD paper
#combines (1) Lark CDL conversion rates (range to range, range to crop) 
# (2) NASS (USDA Census) land rents
# (3) Control variables, including a. road density, b. annual county population c. county farm demographics d. area of different land uses e. irrigated land area
# Land conversion data is for private grassland (not forest, developed or water, not public land)

#Follows on from PerfectEnemyofGood_formatNASSRentdata.r PerfectEnemyofGood_formatUSDAIrrigateddata.r and precursor to PerfectEnemyofGood__econometricmodel.r

#################################################################################################################
#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/Inputs/""
wd <- "Y:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/"
setwd(wd)
options(stringsAsFactors=FALSE) # turn off automatic factor coersion
options(scipen=9999)            # turn off plotting axis lables in scientific notation

#############################
#allDat <- read.csv("Land Ownership and Protected Areas/SGCounties/SGCounties_ownership_habitat_area.csv", header=TRUE) #Only sage grouse counties
sgCounties <- read.csv(paste0(dirname(wd),"/2 Data by county/SGCounty_alldata_160715.csv")) #I only use this dataset to pull the counties that are in sage grouse habitat

#######################
###CONVERSION PROBABILITIES FROM LARK CDL BY LCC
#######################
#Annual conversion probabilites calculated from remotely sensed data on area of crop/range/pasture from 2008-2012
#For grassland (not forest, urban or water) on private land
#data from Tyler Lark (rangeland, cropland)
#Land Capability Class from USDA

#area of cropland from LARK et al. 2015
	croplandAreaDF <- read.csv("Cropland and pasture/Cropland and Rangeland by year/Cropland_Area_byYear_ha_2008_2015.csv")
	croplandAreaAllcounties <- croplandAreaDF[, c("ADMIN_FIPS", "STATE", "NAME", "Area_ha", "X2008_Cropland_area_ha", "X2009_Cropland_area_ha", "X2010_Cropland_area_ha", "X2011_Cropland_area_ha", "X2012_Cropland_area_ha", "X2013_Cropland_area_ha", "X2014_Cropland_area_ha", "X2015_Cropland_area_ha", "X2009_RangetoCropland_area_ha", "X2010_RangetoCropland_area_ha", "X2011_RangetoCropland_area_ha", "X2012_RangetoCropland_area_ha", "X2013_RangetoCropland_area_ha", "X2014_RangetoCropland_area_ha", "X2015_RangetoCropland_area_ha")] #select only relevant columns
	croplandAreaAllcounties <- croplandAreaAllcounties[order(croplandAreaAllcounties$ADMIN_FIPS),] #order data by fips

#area of rangeland from LARK et al. 2015
	rangelandAreaDF <- read.csv("Cropland and pasture/Cropland and Rangeland by year/Rangeland_Area_byYear_ha_2008to2015.csv")
	rangelandAreaAllcounties <- rangelandAreaDF[ , c("ADMIN_FIPS", "Stable.area.of.noncrop", "Rangeland_2008_ha", "Rangeland_2009_ha", "Rangeland_2010_ha", "Rangeland_2011_ha", "Rangeland_2012_ha", "Rangeland_2013_ha", "Rangeland_2014_ha", "Rangeland_2015_ha", "X2009_RangetoRange_ha", "X2010_RangetoRange_ha", "X2011_RangetoRange_ha", "X2012_RangetoRange_ha", "X2013_RangetoRange_ha", "X2014_RangetoRange_ha", "X2015_RangetoRange_ha", "TotalArea_RangeorCrop_ha")] #select only relevant columns
	rangelandAreaAllcounties <- rangelandAreaAllcounties[order(rangelandAreaAllcounties$ADMIN_FIPS),] #order data by fips
CropRangeDF <- merge(croplandAreaAllcounties, rangelandAreaAllcounties, by="ADMIN_FIPS", all.x=TRUE)


###BY YEAR #FOR LCC1to4, 1to6 vs 5or6
LCC1to4RangeDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_LCC1to4.csv", header=TRUE)
LCC1to4CropDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_LCC1to4.csv", header=TRUE)
LCC1to4RangeDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_2015_LCC1to4.csv", header=TRUE)
LCC1to4CropDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_2015_LCC1to4.csv", header=TRUE)
LCC1to4Merge <- merge(LCC1to4RangeDF, LCC1to4CropDF[, !names(LCC1to4CropDF) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC1to4Mergenew <- merge(LCC1to4RangeDFnew, LCC1to4CropDFnew[, !names(LCC1to4CropDFnew) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC1to4Merge <- merge(LCC1to4Merge, LCC1to4Mergenew[, !names(LCC1to4Mergenew) %in% c("Total_Area_RangeorCrop_inLCC1to4", "Stable_cropland")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)

LCC1to6RangeDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_LCC1to6.csv", header=TRUE)
LCC1to6CropDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_LCC1to6.csv", header=TRUE)
LCC1to6RangeDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_2015_LCC1to6.csv", header=TRUE)
LCC1to6CropDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_2015_LCC1to6.csv", header=TRUE)
LCC1to6Merge <- merge(LCC1to6RangeDF, LCC1to6CropDF[, !names(LCC1to6CropDF) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC1to6Mergenew <- merge(LCC1to6RangeDFnew, LCC1to6CropDFnew[, !names(LCC1to6CropDFnew) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC1to6Merge <- merge(LCC1to6Merge, LCC1to6Mergenew[, !names(LCC1to6Mergenew) %in% c("Total_Area_RangeorCrop_inLCC1to6", "Stable_cropland", "Stable.area.of.noncrop")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)

LCC5or6RangeDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_LCC5or6.csv", header=TRUE)
LCC5or6CropDF <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_LCC5or6.csv", header=TRUE)
LCC5or6RangeDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/Rangeland_Area_byYear_2015_LCC5or6.csv", header=TRUE)
LCC5or6CropDFnew <- read.csv("Cropland and pasture/Cropland and Rangeland by LCC/CroplandArea_byYear_2015_LCC5or6.csv", header=TRUE)
LCC5or6Merge <- merge(LCC5or6RangeDF, LCC5or6CropDF[, !names(LCC5or6CropDF) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC5or6Mergenew <- merge(LCC5or6RangeDFnew, LCC5or6CropDFnew[, !names(LCC5or6CropDFnew) %in% c("NAME")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)
LCC5or6Merge <- merge(LCC5or6Merge, LCC5or6Mergenew[, !names(LCC5or6Mergenew) %in% c("Total_Area_RangeorCrop_inLCC5or6", "Stable_cropland", "Stable.area.of.noncrop")], by=c("ADMIN_FIPS", "STATE"), all.x=TRUE)

CropRangeLCCLong <- rbind(with(CropRangeDF,	
								data.frame(
									ADMIN_FIPS=rep(ADMIN_FIPS, 7), 
									State=rep(STATE, 7), 
									Year=rep(c(2009, 2010, 2011, 2012, 2013, 2014, 2015), each=nrow(CropRangeDF)),	
									TotalRangeorCropAreainLCC_ha = rep(TotalArea_RangeorCrop_ha, 7),
									LCC = rep("AllLCC", 7*nrow(CropRangeDF)),
									StableArea_noncrop= rep(Stable.area.of.noncrop, 7),
									Area_croplandt0 = c(X2008_Cropland_area_ha, X2009_Cropland_area_ha, X2010_Cropland_area_ha, X2011_Cropland_area_ha, X2012_Cropland_area_ha, X2013_Cropland_area_ha, X2014_Cropland_area_ha),
									Area_Ranget0_to_Cropt1= c(X2009_RangetoCropland_area_ha, X2010_RangetoCropland_area_ha, X2011_RangetoCropland_area_ha, X2012_RangetoCropland_area_ha, X2013_RangetoCropland_area_ha, X2014_RangetoCropland_area_ha, X2015_RangetoCropland_area_ha),
									Area_Ranget0_to_Ranget1= c(X2009_RangetoRange_ha, X2010_RangetoRange_ha, X2011_RangetoRange_ha, X2012_RangetoRange_ha, X2013_RangetoRange_ha, X2014_RangetoRange_ha, X2015_RangetoRange_ha), 
									Area_Ranget0 = c(Rangeland_2008_ha, Rangeland_2009_ha, Rangeland_2010_ha, Rangeland_2011_ha, Rangeland_2012_ha, Rangeland_2013_ha, Rangeland_2014_ha)
									)),
							with(LCC1to6Merge,	
								data.frame(
									ADMIN_FIPS=rep(ADMIN_FIPS, 7), 
									State=rep(STATE, 7), 
									Year=rep(c(2009, 2010, 2011, 2012, 2013, 2014, 2015), each=nrow(CropRangeDF)),	
									TotalRangeorCropAreainLCC_ha = rep(Total_Area_RangeorCrop_inLCC1to6, 7),
									LCC = rep("LCC1to6", 7*nrow(LCC1to6Merge)),
									StableArea_noncrop= rep(Stable.area.of.noncrop, 7),
									Area_croplandt0 = c(X2008_Cropland_area_ha, X2009_Cropland_area_ha, X2010_Cropland_area_ha, X2011_Cropland_area_ha, X2012_Cropland_area_ha, X2013_Cropland_area_ha, X2014_Cropland_area_ha),
									Area_Ranget0_to_Cropt1= c(X2009_RangetoCropland_area_ha, X2010_RangetoCropland_area_ha, X2011_RangetoCropland_area_ha, X2012_RangetoCropland_area_ha, X2013_RangetoCropland_area_ha, X2014_RangetoCropland_area_ha, X2015_RangetoCropland_area_ha),
									Area_Ranget0_to_Ranget1= c(X2009_RangetoRange_ha, X2010_RangetoRange_ha, X2011_RangetoRange_ha, X2012_RangetoRange_ha, X2013_RangetoRange_ha, X2014_RangetoRange_ha, X2015_RangetoRange_ha),
									Area_Ranget0 = c(Rangeland_2008_ha, Rangeland_2009_ha, Rangeland_2010_ha, Rangeland_2011_ha, Rangeland_2012_ha.y, Rangeland_2013_ha, Rangeland_2014_ha)
									)),
							with(LCC5or6Merge, 
								data.frame(
									ADMIN_FIPS=rep(ADMIN_FIPS, 7), 
									State=rep(STATE, 7), 
									Year=rep(c(2009, 2010, 2011, 2012, 2013, 2014, 2015), each=nrow(CropRangeDF)),	
									TotalRangeorCropAreainLCC_ha = rep(Total_Area_RangeorCrop_inLCC5or6, 7),
									LCC = rep("LCC5or6", 7*nrow(LCC5or6Merge)),
									StableArea_noncrop= rep(Stable.area.of.noncrop, 7),
									Area_croplandt0 = c(X2008_Cropland_area_ha, X2009_Cropland_area_ha, X2010_Cropland_area_ha, X2011_Cropland_area_ha, X2012_Cropland_area_ha, X2013_Cropland_area_ha, X2014_Cropland_area_ha),
									Area_Ranget0_to_Cropt1= c(X2009_RangetoCropland_area_ha, X2010_RangetoCropland_area_ha, X2011_RangetoCropland_area_ha, X2012_RangetoCropland_area_ha, X2013_RangetoCropland_area_ha, X2014_RangetoCropland_area_ha, X2015_RangetoCropland_area_ha),
									Area_Ranget0_to_Ranget1= c(X2009_RangetoRange_ha, X2010_RangetoRange_ha, X2011_RangetoRange_ha, X2012_RangetoRange_ha, X2013_RangetoRange_ha, X2014_RangetoRange_ha, X2015_RangetoRange_ha),
									Area_Ranget0 = c(Rangeland_2008_ha, Rangeland_2009_ha, Rangeland_2010_ha, Rangeland_2011_ha, Rangeland_2012_ha.y, Rangeland_2013_ha, Rangeland_2014_ha)
									)),						
							with(LCC1to4Merge, 
								data.frame(
									ADMIN_FIPS=rep(ADMIN_FIPS, 7), 
									State=rep(STATE, 7), 
									Year=rep(c(2009, 2010, 2011, 2012, 2013, 2014, 2015), each=nrow(CropRangeDF)),	
									TotalRangeorCropAreainLCC_ha = rep(Total_Area_RangeorCrop_inLCC1to4, 7),
									LCC = rep("LCC1to4", 7*nrow(LCC1to4Merge)),
									StableArea_noncrop= rep(Stable_noncrop, 7),
									Area_croplandt0 = c(X2008_Cropland_area_ha, X2009_Cropland_area_ha, X2010_Cropland_area_ha, X2011_Cropland_area_ha, X2012_Cropland_area_ha, X2013_Cropland_area_ha, X2014_Cropland_area_ha),
									Area_Ranget0_to_Cropt1= c(X2009_RangetoCropland_area_ha, X2010_RangetoCropland_area_ha, X2011_RangetoCropland_area_ha, X2012_RangetoCropland_area_ha, X2013_RangetoCropland_area_ha, X2014_RangetoCropland_area_ha, X2015_RangetoCropland_area_ha),
									Area_Ranget0_to_Ranget1= c(X2009_RangetoRange_ha, X2010_RangetoRange_ha, X2011_RangetoRange_ha, X2012_RangetoRange_ha, X2013_RangetoRange_ha, X2014_RangetoRange_ha, X2015_RangetoRange_ha),
									Area_Ranget0 = c(Rangeland_2008_ha, Rangeland_2009_ha, Rangeland_2010_ha, Rangeland_2011_ha, Rangeland_2012_ha.y, Rangeland_2013_ha, Rangeland_2014_ha)
									))
									
								)
CropRangeLCCLong$PercentArea_CropStatic <- with(CropRangeLCCLong, 1-(StableArea_noncrop/TotalRangeorCropAreainLCC_ha))
CropRangeLCCLong$PercentArea_Crop <- with(CropRangeLCCLong, Area_croplandt0/TotalRangeorCropAreainLCC_ha)	
CropRangeLCCLong[CropRangeLCCLong$TotalRangeorCropAreainLCC_ha==0, "PercentArea_Crop"] <- NA #NA where there is no land in that LCC category in that county

##########################
#Calculate and transform RESPONSE VARIABLES
##########################
								
CropRangeLCCLong$ConversionPropCropRangeV2 <- with(CropRangeLCCLong, Area_Ranget0_to_Cropt1/Area_Ranget0_to_Ranget1)
CropRangeLCCLong$ConversionPropCropRangeV2logZIF1 <- log(CropRangeLCCLong$ConversionPropCropRangeV2 + 1)
CropRangeLCCLong$ConversionPropCropRangeV2logZIF001 <- log(CropRangeLCCLong$ConversionPropCropRangeV2 + 0.001)
CropRangeLCCLong$ConversionPropCropRangeV2logZIF0001 <- log(CropRangeLCCLong$ConversionPropCropRangeV2 + 0.0001)
smallval <- min(CropRangeLCCLong[which(CropRangeLCCLong$ConversionPropCropRangeV2>0 & CropRangeLCCLong$LCC=="LCC1to6"), "ConversionPropCropRangeV2"])/2 #half smallest non-zero value
CropRangeLCCLong$ConversionPropCropRangeV2logZIFsmall <- log(CropRangeLCCLong$ConversionPropCropRangeV2 + smallval)

#http://robjhyndman.com/hyndsight/transformations/

write.csv(CropRangeLCCLong, "Cropland and pasture/Allcounties_LarkConversion_byLCC_GrasslandPrivateArea.csv", row.names=FALSE)

#######################
#Merge with NASS RENTAL RATES 
#######################
#created using #PerfectEnemyofGood_formatNASSRentdata.r

rentalRatesNASS <- read.csv("Rental Rates/USDA NASS land rent by county/NASS_LandRents_2008_2014_allRents_long_CPIadjusted.csv", header=TRUE) 
names(rentalRatesNASS) <- c("State_FIPS", "District", "County", "Name", "NASS_FIPS", "Yeart", "NASS_LandRent_Rangeland", "NASS_LandRent_NonIrrigatedCropland", "NASS_DeltaRentCropRangeBetweenYears_t0_t1", "NASS_DeltaRentCropRange")
rentalRatesNASS$NASS_LandRent_RangelandNEG <- rentalRatesNASS$NASS_LandRent_Rangeland*-1
CropRangeLCCLong <- merge(CropRangeLCCLong, rentalRatesNASS, by.y=c("NASS_FIPS", "Yeart"), by.x=c("ADMIN_FIPS", "Year"), all.x=TRUE)

#define counties which overlap sage grouse habitat
CropRangeLCCLong$SageGrouseCounty <- 0
CropRangeLCCLong[CropRangeLCCLong$ADMIN_FIPS %in% sgCounties$ADMIN_FIPS, "SageGrouseCounty"] <- 1
write.csv(CropRangeLCCLong, paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted.csv"), row.names=FALSE)


#########################
#Merge with CONTROL VARIABLES
#########################
CropRangeLCCLong <- read.csv(paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted.csv"), header=TRUE)
controlVariablesDF <- CropRangeLCCLong

###TIME-VARYING
# #2. Add population change - between 2010 & 2012
# library(tidyr)
# popnDFto2009 <- read.csv("Population_USCensus/co-est2009-alldata.csv", header=TRUE, stringsAsFactors=FALSE)
# popnDFto2009$ADMIN_FIPS <- paste0(popnDFto2009$STATE, sapply(popnDFto2009$COUNTY, function(x) formatC(x, width = 3, format = "d", flag = "0")))
# popnDFto2016 <- read.csv("Population_USCensus/co-est2016-alldata.csv", header=TRUE, stringsAsFactors=FALSE)
# popnDFto2016$ADMIN_FIPS <- paste0(popnDFto2016$STATE, sapply(popnDFto2016$COUNTY, function(x) formatC(x, width = 3, format = "d", flag = "0")))
# popnDF <- merge(popnDFto2009, popnDFto2016, by="ADMIN_FIPS", all=TRUE) 
# popnDF1 <- popnDF[,c("ADMIN_FIPS", paste0("POPESTIMATE", 2008:2014))]
# popnDF2 <- popnDF[, c("ADMIN_FIPS", paste0("NPOPCHG_", 2008:2014))]
# names(popnDF1) <- c("ADMIN_FIPS", 2008:2014)
# names(popnDF2) <- c("ADMIN_FIPS", 2008:2014)
# popn1 <- gather(popnDF1, "Year", "Popn", 2:8)
# popn2 <- gather(popnDF2, "Year", "Popn_Chg", 2:8)
# popnDF <- merge(popn1, popn2, by=c("ADMIN_FIPS", "Year"), all=TRUE)	
# popnDF$PopnChg_Perc <- 100*popnDF$Popn_Chg/popnDF$Popn
# popnDF$nextYear = as.numeric(popnDF$Year)+1
# names(popnDF)[2] <- "censusYear"
# write.csv(popnDF, "Population_USCensus/Population_by_county_2008to2015.csv")	
# controlVariablesDF <- merge(controlVariablesDF, popnDF, by.x=c("ADMIN_FIPS","Year"), by.y=c("ADMIN_FIPS","nextYear"), all.x=TRUE)

#STATIC
library(tidyr)
popnDFto2009 <- read.csv("Population_USCensus/co-est2009-alldata.csv", header=TRUE, stringsAsFactors=FALSE)
popnDFto2009$ADMIN_FIPS <- paste0(popnDFto2009$STATE, sapply(popnDFto2009$COUNTY, function(x) formatC(x, width = 3, format = "d", flag = "0")))
popnDFto2016 <- read.csv("Population_USCensus/co-est2016-alldata.csv", header=TRUE, stringsAsFactors=FALSE)
popnDFto2016$ADMIN_FIPS <- paste0(popnDFto2016$STATE, sapply(popnDFto2016$COUNTY, function(x) formatC(x, width = 3, format = "d", flag = "0")))
popnDF <- merge(popnDFto2009, popnDFto2016, by="ADMIN_FIPS", all=TRUE) 
popnDF1 <- popnDF[,c("ADMIN_FIPS", paste0("POPESTIMATE", 2008:2014))]
names(popnDF1) <- c("ADMIN_FIPS", 2008:2014)
popnDF1$Popn <- apply(popnDF1[,c("2008", "2009", "2010", "2011", "2012")], 1, mean)
popnDF2 <- popnDF[, c("ADMIN_FIPS", paste0("NPOPCHG_", 2008:2014))]
names(popnDF2) <- c("ADMIN_FIPS", 2008:2014)
popnDF2$Popn_Chg <- apply(popnDF2[,c("2008", "2009", "2010", "2011", "2012")], 1, mean)
popnDF3 <- merge(popnDF1, popnDF2, by=c("ADMIN_FIPS"), all=TRUE)
popnDF3$PopnChg_Perc <- 100*popnDF3$Popn_Chg/popnDF3$Popn
write.csv(popnDF3, "Population_USCensus/Population_by_county_2008to2015_static.csv")
controlVariablesDF <- merge(controlVariablesDF, popnDF3[,c("ADMIN_FIPS", "Popn", "Popn_Chg", "PopnChg_Perc")], by=c("ADMIN_FIPS"), all.x=TRUE)

###STATIC
#3. Add area of irrigated (harvested ie cropland, not pastureland) land in each county
# This is only available for whole county, not by lcc
#Note the area irrigated is also available in the agDat, though this is the percent irrigated for all agricultural land (cropland, rangeland & pasture)
irrigatedAreaDF <- read.csv("USDA_IrrigatedArea/CroplandArea_Irrigation_USDACensus_byCounty_formatted.csv", header=TRUE)
controlVariablesDF <- merge(controlVariablesDF, irrigatedAreaDF, by="ADMIN_FIPS",all.x=TRUE)

###STATIC
#4. Add Area & % urban (urbanised plus urban clusters)
#data from US Census http://www.census.gov/geo/reference/ua/urban-rural-2010.html
urbanDF <- read.csv("USCensus_urban/PctUrbanRural_County.csv", header=TRUE)
urbanDF$ADMIN_FIPS <- paste0(urbanDF$STATE, sapply(urbanDF$COUNTY, function(x) formatC(x, width = 3, format = "d", flag = "0")))
urbanDF$AREA_URBAN_ha <- urbanDF$AREA_URBAN/10000
controlVariablesDF <- merge(controlVariablesDF, urbanDF[,c("ADMIN_FIPS","AREA_URBAN_ha","AREAPCT_URBAN")], by="ADMIN_FIPS", all.x=TRUE)

###STATIC
#1. Road density
roadDensityDF <- read.csv("Road density/Road_density_byCounty_privatelandsonly_CRedits.csv", header=TRUE)
controlVariablesDF <- merge(controlVariablesDF, roadDensityDF[,c("ADMIN_FIPS", "road_density")], by="ADMIN_FIPS",all.x=TRUE)

#5. Add Land areas in cropland 2012, and other farm demographics
#USDA/NASS ag census data - summarised using Access - AgCensus2012.accdb
agDat <- read.csv("USDA AgCensus/AllCountyData_AgCensus2012.csv", header=TRUE, stringsAsFactors=FALSE)
agNameLookup <- read.csv("USDA AgCensus/Ag_Census_Map_data_variablelookup_07172015.csv", header=TRUE, stringsAsFactors=FALSE) 
#rename columns
for (i in 2:ncol(agDat)){
	names(agDat)[i] <- agNameLookup[which(agNameLookup$MapID == strsplit(names(agDat)[i], "_valueNumeric")[[1]][1]), "MAPTITLE"]
	} 

controlVariablesDF <- merge(controlVariablesDF, agDat, by.x="ADMIN_FIPS", by.y="FipsNumeric", all.x=TRUE)

###STATIC
#6. Add AUM
AUMdf <- read.csv("AUM stocking rates/forageZonal_11states_grassland_byLCC/AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_final.csv")
AUMdf$AUM_propPublicAUM_inCounty <- AUMdf$AUM_public/AUMdf$AUM_total_privateandpublic
AUMdf[is.na(AUMdf$AUM_private),"AUM_propPublicAUM_inCounty"] <- 1 #one where there is not private land
AUMdf[is.na(AUMdf$AUM_public),"AUM_propPublicAUM_inCounty"] <- 0 #zero where there is not public land

controlVariablesDF <- merge(controlVariablesDF, AUMdf, all.x=TRUE, by=c("ADMIN_FIPS", "LCC") )

######################
#Merge CONTROL VARIABLES with RESPONSE VARIABLES
#####################

write.csv(controlVariablesDF, paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static.csv"), row.names=FALSE)
#write.csv(controlVariablesDF, paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_timevarying.csv"), row.names=FALSE)


######################
#add field describing SGMZ county falls within
######################
controlVariablesDF <- read.csv(paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static.csv"), header=TRUE)
sgmz <- read.csv("Boundaries/SGCounties/SGCounties_overlapping_PACSandBreeding_bySGMZregion.csv", header=TRUE)

sgmzDF <- merge(controlVariablesDF, sgmz[,c("ADMIN_FIPS", "Region")], by="ADMIN_FIPS", all.x=TRUE)
write.csv(sgmzDF, paste0(dirname(wd), "/2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static_plus_sgmzs.csv"), row.names=FALSE)
