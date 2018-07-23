#
library(tidyr) #for spread() which is like reshape

options(stringsAsFactors=FALSE) # turn off automatic factor coersion

wd <- "C:/Claire/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)


#######################
###SET UP DATA
#######################
#countyList <- read.csv("1 Inputs/Boundaries/county_list_r.csv") #not using this - instead pull fips from SGCounty_alldata_160715.csv

#######################
###CONVERSION PROBABILITIES FROM NRI
###LAND RENTS FROM NASS CENSUS
#######################

allDat <- read.csv("2 Data by county/SGCounty_alldata_160715.csv") #the rental rates in this dataset come from NASS and the conversion probabilities from the NRI, and I use alternate sources of rental rate and conversion data in some models below; but it also contains data on land cover, AUM, area of sage grouse habitat, area of private/public land, farm characteristics for each county.

######################
###ECONOMIC - ANNUAL LAND RENT DATA from NASS
#######################
#Annual land rent data from NASS available from 2008 onwards
#RENT, CASH, PASTURELAND - EXPENSE, MEASURED IN $ / ACRE etc

	rentDat <- read.csv(paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/raw/NASS_LandRents_2008_2016.csv"), header=TRUE, stringsAsFactors=TRUE)
		rentDat$County <- as.character(rentDat$County)
	
	Rents <- list(	Pasture_rent="RENT, CASH, PASTURELAND - EXPENSE, MEASURED IN $ / ACRE",
				NonIrrigatedCropland_rent="RENT, CASH, CROPLAND, NON-IRRIGATED - EXPENSE, MEASURED IN $ / ACRE", 
				IrrigatedCropland_rent="RENT, CASH, CROPLAND, IRRIGATED - EXPENSE, MEASURED IN $ / ACRE")
	Years <- c(2008:2014) #there is no data for 2015
	
	#Preprocessing of NASS land rent data
	#function to convert caps to camel case
		camel <- function(x) {
			x <- tolower(x)
			s <- strsplit(x, " ")[[1]]
			paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
		}
	#convert county names in rent dataset to camel case
		rentDat$NewCounty <- sapply(rentDat$County, function(x) camel(x))

	#subset data
		rentDat <- rentDat[rentDat$Data.Item %in% unlist(Rents) & rentDat$Year %in% Years ,c("State", "State.ANSI", "Ag.District", "Ag.District.Code", "NewCounty", "County", "County.ANSI", "Year", "Data.Item", "Value")]
		rentDat$Value <- as.numeric(gsub(",", "", rentDat$Value)) #remove commas from values and convert to numeric

	#Adjust for CPI
		CPItable <- read.csv(paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/CPI_adjustment.csv"), header=TRUE)
		rentDat <- merge(rentDat, CPItable, by="Year")
		rentDat$Value2015 <- round(rentDat$Value*rentDat$Dollars2015, 2)
		rentDat <- rentDat[, !(colnames(rentDat) %in% c("Value", "CPI", "Dollars2012", "Dollars2015"))] #drop columns
		
	###Reshape data
		rentWide <- spread(rentDat, Year, value=Value2015, fill=NA)
		rentWide$Av <- round(rowMeans(rentWide[,dput(as.character(Years))], na.rm=TRUE), 2)
		rentWide$Min <- apply(rentWide[,dput(as.character(Years))],1, min, na.rm=TRUE)
		rentWide$Max <- apply(rentWide[,dput(as.character(Years))],1, max, na.rm=TRUE)
		head(rentWide)
		names(rentWide)[which(names(rentWide)=="State")] <- "STATE"
		names(rentWide)[which(names(rentWide)=="County")] <- "COUNTY"

	#save reshaped data
		write.csv(rentWide, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/intermediate/NASS_LandRents_", min(Years), "_", max(Years), "_reshaped_CPIadjusted.csv"), row.names=FALSE)

	###Fill in gaps for counties with no data
	#Table of which district counties are in to fill in gaps
		countyList <- read.csv(paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/county_list_r.csv"), stringsAsFactors=FALSE)
		countyList <- countyList[countyList$Flag==1,]#pull only current counties
		countyList <- countyList[countyList$County<888 & countyList$County>0,]#exclude combined counties & state totals
		countyList <- countyList[countyList$County<888,]#exclude combined counties & state totals
		countyList$NASS_FIPS <- paste0(countyList$State, formatC(countyList$County, width=3, format="d", flag="0"))
		#Table of data to fill with
		#filldf <- rentWide[rentWide$COUNTY=="OTHER (COMBINED) COUNTIES",]

	#Make a list of dataframes (one list item per Rent)
	rentList <- lapply(Rents, function(currRent){
		currSub <- droplevels(rentWide[rentWide$Data.Item==currRent,])
		currdf <- merge(countyList[,1:4], currSub, by.x=c("Name", "State"), by.y=c("NewCounty", "State.ANSI"), all.x=TRUE)
			
		#fill counties with no data using the average for the Ag.District they fall within
		filldf <- currSub[currSub$COUNTY=="OTHER (COMBINED) COUNTIES",]#Table of which values to fill with
		nodata <- which(is.na(currdf$Av)) #find gaps
		print(names(Rents[which(Rents %in% currRent)]))
		length(nodata)
		#print(nodata)
		
		#print(currdf[nodata,1:3])
		
		#fill in rows
			for(i in nodata){
				currRow <- currdf[i,]
				currFill <- filldf[	filldf$Ag.District.Code==currRow[["District"]] & 
									filldf$State.ANSI==currRow[["State"]] ,
									c("STATE", "Ag.District", "Ag.District.Code", "COUNTY", "County.ANSI", "Data.Item", as.character(Years), "Av", "Min", "Max")	]
				if(nrow(currFill)>0){
				currdf[i,]<- data.frame(currdf[i, 1:4], currFill)
				}
			}
		write.csv(currdf, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/intermediate/NASS_LandRents_", min(Years), "_", max(Years), "_", names(Rents[which(Rents %in% currRent)]), "_CPIadjusted.csv"), row.names=FALSE) 
		return(currdf)
	})
	
	#output dataframe in wide format
		rentDFwide <- data.frame(countyList, do.call(cbind, lapply(rentList, data.frame, stringsAsFactors=FALSE)))
		rentDFwide <- rentDFwide[, c(1:6, 17:26, 37:46, 57:66)]

		
		#Calculate deltarents within years
		rentDFwide$DeltaRentRangeCrop2008 <- with(rentDFwide, NonIrrigatedCropland_rent.X2008 - Pasture_rent.X2008)
		rentDFwide$DeltaRentRangeCrop2009 <- with(rentDFwide, NonIrrigatedCropland_rent.X2009 - Pasture_rent.X2009)
		rentDFwide$DeltaRentRangeCrop2010 <- with(rentDFwide, NonIrrigatedCropland_rent.X2010 - Pasture_rent.X2010)
		rentDFwide$DeltaRentRangeCrop2011 <- with(rentDFwide, NonIrrigatedCropland_rent.X2011 - Pasture_rent.X2011)
		rentDFwide$DeltaRentRangeCrop2012 <- with(rentDFwide, NonIrrigatedCropland_rent.X2012 - Pasture_rent.X2012)
		rentDFwide$DeltaRentRangeCrop2013 <- with(rentDFwide, NonIrrigatedCropland_rent.X2013 - Pasture_rent.X2013)
		rentDFwide$DeltaRentRangeCrop2014 <- with(rentDFwide, NonIrrigatedCropland_rent.X2014 - Pasture_rent.X2014)
		#rentDFwide$DeltaRentRangeCrop2015 <- with(rentDFwide, NonIrrigatedCropland_rent.X2015 - Pasture_rent.X2015)
		#Calculate deltarents across years
		rentDFwide$DeltaRent08to09 <- with(rentDFwide, DeltaRentRangeCrop2009 - DeltaRentRangeCrop2008)
		rentDFwide$DeltaRent09to10 <- with(rentDFwide, DeltaRentRangeCrop2010 - DeltaRentRangeCrop2009)
		rentDFwide$DeltaRent10to11 <- with(rentDFwide, DeltaRentRangeCrop2011 - DeltaRentRangeCrop2010)
		rentDFwide$DeltaRent11to12 <- with(rentDFwide, DeltaRentRangeCrop2012 - DeltaRentRangeCrop2011)
		rentDFwide$DeltaRent12to13 <- with(rentDFwide, DeltaRentRangeCrop2013 - DeltaRentRangeCrop2012)
		rentDFwide$DeltaRent13to14 <- with(rentDFwide, DeltaRentRangeCrop2014 - DeltaRentRangeCrop2013)
		#rentDFwide$DeltaRent14to15 <- with(rentDFwide, DeltaRentRangeCrop2015 - DeltaRentRangeCrop2014)	
			
	write.csv(rentDFwide, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/NASS_LandRents_", min(Years), "_", max(Years), "_allRents_wide_CPIadjusted.csv"), row.names=FALSE)
	
	#output data for only sage grouse counties
		rentDFwideSG <- rentDFwide[rentDFwide$NASS_FIPS %in% allDat$ADMIN_FIPS,]
		rentDFwideSG <- rentDFwideSG[order(rentDFwideSG$NASS_FIPS),] #order by FIPS
	write.csv(rentDFwide, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/NASS_LandRents_", min(Years), "_", max(Years), "_allRents_wide_sgcounties_CPIadjusted.csv"), row.names=FALSE)

	
	#output dataframe in long format
		rentDFlong <- with(rentDFwide, data.frame(State=rep(State, 7), District=rep(District, 7), County=rep(County, 7), Name=rep(Name, 7), NASS_FIPS=rep(NASS_FIPS, 7), Yeart=rep(Years, each=nrow(rentDFwide)), LandRent_Rangeland=unlist(rentDFwide[, c("Pasture_rent.X2008", "Pasture_rent.X2009", "Pasture_rent.X2010", "Pasture_rent.X2011", "Pasture_rent.X2012","Pasture_rent.X2013","Pasture_rent.X2014")], use.names=FALSE), LandRent_NonIrrigatedCropland=unlist(rentDFwide[,c("NonIrrigatedCropland_rent.X2008", "NonIrrigatedCropland_rent.X2009", "NonIrrigatedCropland_rent.X2010", "NonIrrigatedCropland_rent.X2011", "NonIrrigatedCropland_rent.X2012", "NonIrrigatedCropland_rent.X2013", "NonIrrigatedCropland_rent.X2014")], use.names=FALSE), DeltaRentCropRangeBetweenYears_t0_t1= c(rep(NA, nrow(rentDFwide)), unlist(rentDFwide[,c("DeltaRent08to09", "DeltaRent09to10", "DeltaRent10to11", "DeltaRent11to12", "DeltaRent12to13", "DeltaRent13to14")]))))
		rentDFlong$DeltaRentCropRange <- with(rentDFlong, LandRent_NonIrrigatedCropland - LandRent_Rangeland)
	write.csv(rentDFlong, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/NASS_LandRents_", min(Years), "_", max(Years), "_allRents_long_CPIadjusted.csv"), row.names=FALSE)
		
	#output data for only sage grouse counties
		rentDFlongSG <- rentDFlong[rentDFlong$NASS_FIPS %in% allDat$ADMIN_FIPS,]
		rentDFlongSG <- rentDFlongSG[order(rentDFlongSG$NASS_FIPS),] #order by FIPS
	write.csv(rentDFlongSG, paste0(wd, "1 Inputs/Rental Rates/USDA NASS land rent by county/NASS_LandRents_", min(Years), "_", max(Years), "_allRents_long_sgcounties_CPIadjusted.csv"), row.names=FALSE)
	
	
