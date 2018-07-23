

#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/Inputs/""
#wd <- "C:/Claire/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/"
wd <- "Y:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/"
setwd(wd)

############################################################################################
#Calculate the proportion of the cropland in a county that is irrigated:
AreaIrrigatedDF <- read.csv("USDA_IrrigatedArea/CroplandArea_Irrigation_USDACensus.csv", header=TRUE, stringsAsFactors=FALSE)
names(AreaIrrigatedDF)[c(19,21,23,25)] <- c("Cropland_Acres", "Cropland_Harvested_Acres", "Cropland_Harvested_Irrigated_Acres", "Cropland_PasturedOnly_Acres")

#remove commas from values and convert to numeric
AreaIrrigatedDF$Cropland_Acres <- as.numeric(gsub(",", "", AreaIrrigatedDF$Cropland_Acres))
AreaIrrigatedDF$Cropland_Harvested_Acres <- as.numeric(gsub(",", "", AreaIrrigatedDF$Cropland_Harvested_Acres))
AreaIrrigatedDF$Cropland_Harvested_Irrigated_Acres <- as.numeric(gsub(",", "", AreaIrrigatedDF$Cropland_Harvested_Irrigated_Acres))
AreaIrrigatedDF$Cropland_PasturedOnly_Acres <- as.numeric(gsub(",", "", AreaIrrigatedDF$Cropland_PasturedOnly_Acres))

AreaIrrigatedDF$PercentCroplandthatisIrrigated <- with(AreaIrrigatedDF, round(Cropland_Harvested_Irrigated_Acres/Cropland_Acres, 4))
AreaIrrigatedDF$PercentCroplandthatisIrrigated[is.infinite(AreaIrrigatedDF$PercentCroplandthatisIrrigated)] <- 0
AreaIrrigatedDF$PercentCroplandthatisIrrigated[is.na(AreaIrrigatedDF$PercentCroplandthatisIrrigated)] <- 0


AreaIrrigatedDF$ADMIN_FIPS <- paste0(AreaIrrigatedDF$State.ANSI, sapply(AreaIrrigatedDF$County.ANSI, function(x) formatC(x, width = 3, format = "d", flag = "0")))

AreaIrrigatedDFOut <- AreaIrrigatedDF[, c("ADMIN_FIPS", "Cropland_Acres", "Cropland_Harvested_Acres", "Cropland_Harvested_Irrigated_Acres", "Cropland_PasturedOnly_Acres", "PercentCroplandthatisIrrigated")]

write.csv(AreaIrrigatedDFOut, "USDA_IrrigatedArea/CroplandArea_Irrigation_USDACensus_byCounty_formatted.csv", row.names=FALSE)
