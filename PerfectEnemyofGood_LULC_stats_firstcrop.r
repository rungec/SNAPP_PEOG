#What is the first crop in a county?

library(plyr)
library(rgdal)

wd <- "N:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics"
setwd(wd)

firstcrop <- read.csv("1 Inputs/Cropland and pasture/First_crop_byCounty_2008to2012.txt", header=TRUE) #from tabulate by area first_crop_class.tif with sgCounties as zones
#set names as class
reclassTbl <- read.csv("1 Inputs/LULC/CDL reclassification/CDL_levels_reclass2.csv", header=TRUE)
classids <- as.numeric(sapply(names(firstcrop)[3:length(names(firstcrop))], function(x) strsplit(x, "VALUE_")[[1]][2]))
newnames <- as.character(reclassTbl[reclassTbl$ID %in% classids, "Class_Names"])
names(firstcrop)[3:length(names(firstcrop))] <- newnames

#summarise firstcrop as %
firstcrop$Wheat <- apply(firstcrop[,c("Winter Wheat","Spring Wheat", "Durum Wheat")], 1, sum)
firstcrop$totalarea <- apply(firstcrop[,4:(ncol(firstcrop)-1)], 1, sum)
firstcrop_percent <- data.frame(firstcrop[1:2], round(100*firstcrop[,4:ncol(firstcrop)]/firstcrop$totalarea, 3))

sgCountiesShp <- readOGR("1 Inputs/Boundaries/SGCounties", "SGCounties_overlapping_PACSandBreeding_NAD83alb", integer64='warn.loss')

	#make a shp with all the data
	firstcropShp <- merge(sgCountiesShp, firstcrop_percent, by.x='fips', by.y='ADMIN_FIPS', all.x=TRUE) 

county_mostcommon <- lapply(1:nrow(firstcropShp@data), function(x){
					most_common_breakout_crop = ifelse(firstcropShp@data$totalarea[x]!=100, "NA", names(which.max(firstcropShp@data[x,17:83])))														
					percent_breakoutcrop = max(firstcropShp@data[x,17:83])
					return(data.frame(most_common_breakout_crop, percent_breakoutcrop))
					})
					

county_mostcommonDF <- do.call(rbind, county_mostcommon)

firstcropShp@data <- data.frame(firstcropShp@data, county_mostcommonDF)
	
write.csv(firstcropShp@data, "1 Inputs/Cropland and pasture/First_crop_byCounty_2008to2012_aspercent.csv")
writeOGR(firstcropShp, "1 Inputs/Cropland and pasture", "First_crop_byCounty_2008to2012_aspercent", driver="ESRI Shapefile", overwrite=TRUE)

