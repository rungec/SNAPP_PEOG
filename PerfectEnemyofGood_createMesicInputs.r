##This script summarises mesic habitat in each county, and how much cropland conversion occurs on mesic habitat (held in private land).
#Files used here were processed in ArcGIS - see PerfectEnemyofGood_Processing_notes.md for details
#the output files are used in PerfectEnemyofGood_model_marginaleffects.R


options(stringsAsFactors=FALSE) # turn off automatic factor coersion
options(scipen=9999) #turn off scientific notation

library(plyr)
library(reshape2)

#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
wd <- "N:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)

####################
#SET UP DATA
####################
mesicDF <- read.csv("1 Inputs/Sage grouse/Mesic_habitat_Donnelly/MESIC_private_Lark_diss.csv", header=TRUE)
#drop blank row
mesicDF <- mesicDF[!is.na(mesicDF$ADMIN_FIPS),]
mesicDF$habitat <- "mesic"
mesicDF$habitat[mesicDF$TYPE=="non_mesic"] <- "non_mesic"
mesicDF$convertedtocrop <- 0
mesicDF$convertedtocrop[mesicDF$gridcode==3] <- 1

mesicLong <- aggregate(Area_ha ~ ADMIN_FIPS + habitat + convertedtocrop, data=mesicDF, sum)
mesicWide <- dcast(mesicLong, formula= ADMIN_FIPS ~ habitat + convertedtocrop, value.var="Area_ha", fun.aggregate=sum)
mesicWide$Area_mesic_ha <- with(mesicWide, mesic_0 + mesic_1)
mesicWide$prop_mesic_converted <- with(mesicWide, mesic_1/(mesic_0 + mesic_1))
mesicWide$prop_mesic_converted[is.na(mesicWide$prop_mesic_converted)] <- 0
mesicWide$prop_nonmesic_converted <- with(mesicWide, non_mesic_1/(non_mesic_0 + non_mesic_1))
mesicWide$prop_nonmesic_converted[is.na(mesicWide$prop_nonmesic_converted)] <- 0
mesicWide$prop_conversion_on_mesic <- with(mesicWide, mesic_1/(mesic_1 + non_mesic_1))
mesicWide$prop_conversion_on_mesic[is.na(mesicWide$prop_conversion_on_mesic)] <- 0


write.csv(mesicWide, "1 Inputs/Sage grouse/Mesic_habitat_Donnelly/MESIC_proportionConverted_byCounty.csv", row.names=FALSE)
