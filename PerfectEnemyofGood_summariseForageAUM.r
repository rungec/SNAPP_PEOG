
library(parallel)  
no_cores <-4
# test=function() {source(paste0(wd, "PerfectEnemyofGood_createInputs.r")}
# test() #to run

options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off plotting axis lables in scientific notation

library(rgdal)
library(rgeos) #for gIsValid
library(plyr)
library(raster)
library(tidyr) #for spread
library(foreign) #for read.dbf

#library(foreach)

########################
## Setting system preferences for better work with raster
# Define your temp folder
# my_tmpdir=paste0(wd, "tmpRaster")

# # Create it (handles the case where the folder already exists)
# dir.create(my_tmpdir, showWarnings=F)

# # Set the raster option to this folder
# rasterOptions(tmpdir= my_tmpdir)

## Maximum Memory
# settings:
# rasterOptions(maxmemory =1e+09)
# rasterOptions(chunksize=1e+08)

#############################
#Set dirs
#wd <- "C:/Claire/"
wd <- "/home/runge/Data/"
US_STATES = c("CO", "ID", "MT" ,"ND" ,"NV", "OR", "SD" ,"UT", "WA" ,"WY", "CA")

########################
# #SET UP SHAPEFILE FOR sg counties, all land, by sage grouse habitat & ownership
# setwd(paste0(wd, "NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/"))
# countyshp <- readOGR("Boundaries/SGCounties/" , "SGCounties_overlapping_PACSandBreeding_NAD83alb")
# dat <- readOGR("Land Ownership and Protected Areas/SGCounties", "SGCounties_sgPACsBr_and_ownership_union_eliminate")
# datRast <- raster("Land Ownership and Protected Areas/SGCounties/SGCounties_sgPACsBr_rast.tif")
# #crs(datRast) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #Check proj is EPGS:5070 Albers Conic Equal Area 
# #for some reason the crs inrasters created in ArcGIS or arcpy is different than the crs of those created in R. The rasters overlap, and still overlap after the projection is assigned in this way in R. This was the best work around we could come up with.
# privatelandownerslist <- c("Private Land", "Private Conservation Land")
# #dput(unique(landowner@data$own_type))
# #Test for geometry problems and try to fix
# fixshp <- function(x){
		# if(gIsValid(x)==FALSE){ # any bad polys?
		# print("Invalid geometry - performing zero buffer")
			# x <- gBuffer(x, byid=TRUE, width=0) #fix
		# } else {
			# x <- x
		# }
# }

# dat <- fixshp(dat)
# countyshp <- fixshp(countyshp)

#if it throws errors
#gIsValid() to check

###########################

#extract gSSURGO data for each polygon in dat
#forage productivity is an estimated value for each soil class and is drawn from the component table within gssurgo
#uses data extracted using scripts JoinMosaic_SGMZ.R & JoinMosaic_SGMZ_gdbExtractor.py
#these scripts are stored in the soil-organic-carbon repository

#Summarise forage for the sage grouse counties, polygons split by owner and sage grouse habitat
#uses a raster (datRast), that is a rasterized version of the zonal shp
# zonalFun <- function(currState){
		# print(currState)
		# dir.create(my_tmpdir, showWarnings=F)
		# rasterOptions(tmpdir= my_tmpdir)
		
		# rastList <- list.files("AUM stocking rates/forage/", pattern=currState, full.names=TRUE)
		# rastStack <- stack(rastList, RAT=FALSE)
		# print(paste("starting crop", currState, as.character(Sys.time()), sep=" "))		
		# datRastSub <- crop(datRast, extent(rastStack[[1]]))
		# rastStack <- crop(rastStack, extent(datRastSub))
		# print(paste("starting zonal", currState, as.character(Sys.time()), sep=" "))
		# forageZonal <- zonal(rastStack, datRastSub, 'mean', na.rm=TRUE)
		# colnames(forageZonal)<- c("FID", sapply(rastList, function(i) paste(strsplit(basename(i), "_")[[1]][2:3], collapse="_")))
		# write.csv(forageZonal, filename=paste0("AUM stocking rates/forageZonalSGCounties_AllLand/", currState, ".csv"))
		# return(forageZonal)
		# print(paste("finished zonal ", currState, as.character(Sys.time()), sep=" "))
		# unlink(my_tmpdir, recursive = TRUE)
		# }

# datforage <- foreach(mystate=US_STATES) %dopar% {zonalFun(mystate)}


##########################
#CALCULATE FORAGE FOR 11states, by LCC & county, grasslands only
##########################
setwd(paste0(wd, "/NCEAS_Postdoc/"))

###set up rasters
lcc <- raster("P1 Sage Grouse/Data/Original/LULC/Land Capability Class/LCC_100m.tif") #land capability class
lulc <- raster("P4 ACR revised methods/Analysis/rasters/CDL raster masked by grasslands/LarkCDL_grassland.tif")

#resample lcc.tif to match lulc
lccresample <- crop(lcc, extent(lulc), snap='near')
lccresample <- resample(lccresample, lulc, method='ngb')
lccresample <- crop(lccresample, extent(lulc))

#reclass lulc.tif to NA if forest, developed or water
lulcreclass <- subs(x=lulc, y=data.frame(by=c(1,2,3,4,5), which=c(1,1,1,1,1)), subswithNA=TRUE)

#combine lulc.tif and lcc.tif
maskLCCLULC <- overlay(lccresample, lulcreclass, fun=function(x,y) { x*y*1000000 })
writeRaster(maskLCCLULC, "P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/LULC/LCC_grasslands_mask_56m.tif", format='GTiff', datatype="INT4S") #the value in this rasters is the land capability class (100000=LCC1, 2000000=LCC2, etc). Value = NA where forest, developed or water

###Summarise forage by county and ownership (public/private) for the 11 states, for land that is not forest, developed or water, and excluding strictly protected areas, by LCC
#each value (integer) in maskRast corresponds to the FID on PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp
setwd(paste0(wd, "/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/1 Inputs/"))

#load combined lcc and lulc raster created above
maskLCCLULC <- raster("/LULC/LCC_grasslands_mask_56m.tif")

#calculate mean forage for each county, split by lcc and public/private
datforage2 <- mclapply(US_STATES, function(currState){
		print(currState)
		
		rastList <- list.files("AUM stocking rates/forage/", pattern=currState, full.names=TRUE)
		rastStack <- stack(rastList, RAT=FALSE)
		
		print(paste("starting create mask", currState, as.character(Sys.time()), sep=" "))
		#create mask
			maskTemplate <- rastStack
			maskTemplate[] <- 0
			writeRaster(maskTemplate, filename=sprintf("AUM stocking rates/intermediate/%s.tif", currState), format="GTiff", datatype="INT4S", overwrite=TRUE)
			template_extent <- paste(maskTemplate@extent@xmin, maskTemplate@extent@ymin, maskTemplate@extent@xmax, maskTemplate@extent@ymax, sep=" ")
			rm(maskTemplate)
		#burn in protected areas - value = FID of PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp
			print(paste0("Starting rasterize 11states", Sys.time()))
			system(paste(sprintf("gdal_rasterize -at -a Id -of GTiff -tr 30 30 -te %s", template_extent),shQuote("Land Ownership and Protected Areas/11states/PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp"), shQuote(sprintf("AUM stocking rates/intermediate/%s.tif", currState))))
		#merge with fid with lcc
			maskPA <- raster(sprintf("AUM stocking rates/intermediate/%s.tif", currState))
			print("startingcrop")
			maskLCC <- crop(maskLCCLULC, extent(maskPA), snap='out')
			maskLCC <- resample(maskLCC, maskPA, method='ngb', filename=sprintf("AUM stocking rates/intermediate/%s_resamplelcc.tif", currState), format="GTiff", datatype="INT4S", overwrite=TRUE)
			print("starting overlay")			
			maskRast <- overlay(maskPA, maskLCC, fun=function(x,y) {x+y})
			writeRaster(maskRast, filename=sprintf("AUM stocking rates/intermediate/%s_mask.tif", currState), format="GTiff", datatype="INT4S", overwrite=TRUE)
			rm(maskLCC)
			rm(maskPA)
		#summarise mean forage for each admin_fips/public.private/lcc zone
		print(paste("starting zonal", currState, as.character(Sys.time()), sep=" "))
			forageZonal <- zonal(rastStack, maskRast, 'mean', na.rm=TRUE)
			write.csv(forageZonal, file=paste0(getwd(), "/AUM stocking rates/forageZonal_11states_grassland_byLCC/", currState, ".csv"), row.names=FALSE)
		return(forageZonal)
		print(paste("finished zonal ", currState, as.character(Sys.time()), sep=" "))
		#unlink(my_tmpdir, recursive = TRUE)
		}, mc.cores=no_cores)

#first digit = LCC; second digit=public(0), private(1); final 5 digits=ADMIN_FIPS

###########################
#calculate area of land in each LCC, how much falls in private or public, by County
###########################
areaDF <- mclapply(US_STATES, function(currState){
	currRast <- raster(sprintf("AUM stocking rates/intermediate/%s_mask.tif", currState))
	areaTbl <- data.frame(freq(currRast))
	areaTbl$State <- rep(currState, nrow(areaTbl))
	return(areaTbl)
	}, mc.cores=no_cores)

areaDFallStates <- do.call(rbind, areaDF)
tabulatedAreas <- ddply (areaDFallStates, c("value", "State"), summarise, Freq_total=sum(count)) 
tabulatedAreas$Area_ha <- round(tabulatedAreas$Freq_total *30*30/10000, 5) #raster is 30m resolution
tabulatedAreas$Area_acres <- round(tabulatedAreas$Area_ha*2.47105, 3)

#The rasters for each state overlapped the edges of neighbouring states. This dataset includes rows that are assigned to these states but where the counties are actually in neighbouring states
write.csv(tabulatedAreas, "LULC/Total_area_byLCC_and_publicprivate_11states.csv", row.names=FALSE)

###########################
#Combine forage data for all states
###########################
forageDir <- "AUM stocking rates/forageZonal_11states_grassland_byLCC/"
#FIDS in forageDFs are drawn from PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp
#first digit = LCC; second digit=public(0), private(1); final 5 digits=ADMIN_FIPS


forageDFList <- lapply(US_STATES, function(currState){
			forageDF <- read.csv(paste0(forageDir, currState, ".csv")) #load forage data for that state
			names(forageDF)[2:4] <- c("rsprod_h", "rsprod_l", "rsprod_r")
			#extract lcc from zone
			forageDF$LCC <- as.numeric(sapply(forageDF$zone, function(x) strsplit(formatC(x, width=7, format="d", flag=0), "")[[1]][1]))
			#extract private or public from zone (1 if private, 0 if public)
			forageDF$Private <- as.numeric(sapply(forageDF$zone, function(x) strsplit(formatC(x, width=7, format="d", flag=0), "")[[1]][2]))
			#extract FID from zone
			forageDF$FID <- forageDF$zone - 1000000*forageDF$LCC
			#extract ADMIN_FIPS from zone
			forageDF$ADMIN_FIPS <- as.numeric(sapply(forageDF$zone, function(x) paste0(strsplit(formatC(x, width=7, format="d", flag=0), "")[[1]][3:7], collapse="")))
			forageDF$State <- rep(currState, nrow(forageDF))
			return(forageDF)
})
combinedforageDF <- do.call(rbind, forageDFList)

###########################
#Calculate AUM
###########################
#Calculate per acre, but not acreage sums as there are NAs in the dataset - which will cause the carrying capacity to be lower than actual - instead, once data is summarised by county then use the county mean to calculate carrying capacity.
AUM <- lapply(list("rsprod_h", "rsprod_r", "rsprod_l"), function(x){
			AUM_peracre = round(combinedforageDF[,x]*0.25/(30*30.5), 3)  #Carrying capacity in AUM #grazin efficiency 25% #30lb air dried feed per day, 30.5days per month
			newcols <- data.frame(AUM_peracre)
			names(newcols) <- paste(names(newcols), strsplit(x, "_")[[1]][2], sep="_")
			return(newcols)
})

combinedforageDF  <- data.frame(combinedforageDF, AUM[[1]], AUM[[2]], AUM[[3]])
write.csv(combinedforageDF, file=paste0(forageDir, "AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_raw.csv"), row.names=FALSE)

##########################
#Summarise AUM for each county
##########################
#drop rows with NAs
combinedforageDF2 <- combinedforageDF [complete.cases(combinedforageDF$rsprod_r), ]
#drop rows with LCC=0
combinedforageDF2 <- combinedforageDF2[combinedforageDF2$LCC != 0, ]

combinedforageDF2$LCC1to6 <- 0
combinedforageDF2[combinedforageDF2$LCC %in% c(1,2,3,4,5,6), "LCC1to6"] <- 1
combinedforageDF2$LCC1to4 <- 0
combinedforageDF2[combinedforageDF2$LCC %in% c(1,2,3,4), "LCC1to4"] <- 1
	 
###
#add in the areas for each fid (calculated above)
tabulatedAreas <- read.csv("LULC/Total_area_byLCC_and_publicprivate_11states.csv")
combinedforageDF2 <- base::merge(combinedforageDF2, tabulatedAreas, all.x=TRUE, by.x=c("zone", "State"), by.y=c("value", "State"))

#Calculate the AUM for each zone (i.e. each combination of county-lcc-private.public)
combinedforageDF2$AUM = round(combinedforageDF2$AUM_peracre_r*combinedforageDF2$Area_acres, 2) 
write.csv(combinedforageDF2, file=paste0(forageDir, "AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_raw_allrows.csv"), row.names=FALSE)

#The rasters for each state overlapped the edges of neighbouring states. Drop the rows that are assigned to these states but where the counties are actually in neighbouring states
countylist <- read.dbf("Boundaries/SGCounties/SGCounties_overlapping_PACSandBreeding_NAD83alb.dbf")
combinedforageDF3 <- base::merge(countylist, combinedforageDF2, by.x=c("fips", "STATE"), by.y=c("ADMIN_FIPS", "State"), all.x=TRUE)

#Summarise by LCC groups
AllLCCSum <- ddply(combinedforageDF3, c("ADMIN_FIPS", "Private"), summarise, AUM_total=sum(AUM))
AllLCCSumWide <- spread(AllLCCSum, Private, AUM_total)
LCC1to6Sum <- ddply(combinedforageDF3, c("ADMIN_FIPS", "Private", "LCC1to6"), summarise, AUM_total=sum(AUM))
LCC1to6SumWide <- spread(LCC1to6Sum, Private, AUM_total)
LCC1to4Sum <- ddply(combinedforageDF3, c("ADMIN_FIPS", "Private", "LCC1to4"), summarise, AUM_total=sum(AUM))
LCC1to4SumWide <- spread(LCC1to4Sum, Private, AUM_total)

#reformat the dataset
finalAUMDF <- rbind(AllLCCSumWide[,c("ADMIN_FIPS", "0", "1")], LCC1to4SumWide[LCC1to4SumWide$LCC1to4==1,c("ADMIN_FIPS", "0", "1")], LCC1to6SumWide[LCC1to6SumWide$LCC1to6==1,c("ADMIN_FIPS", "0", "1")], LCC1to6SumWide[LCC1to6SumWide$LCC1to6==0,c("ADMIN_FIPS", "0", "1")])
finalAUMDF$LCC <- c(rep("AllLCC", nrow(AllLCCSumWide)), rep("LCC1to4", nrow(LCC1to4SumWide[LCC1to4SumWide$LCC1to4==1,])), rep("LCC1to6", nrow(LCC1to6SumWide[LCC1to6SumWide$LCC1to6==1,])), rep("LCC7or8", nrow(LCC1to6SumWide[LCC1to6SumWide$LCC1to6==0,])))
names(finalAUMDF)[2:3] <- c("AUM_private", "AUM_public")
finalAUMDF <- finalAUMDF[,c("ADMIN_FIPS", "LCC", "AUM_private", "AUM_public")]

#calculate ratio of public AUM to private AUM
finalAUMDF$AUM_ratio_privatetopublic <- finalAUMDF$AUM_private/finalAUMDF$AUM_public
finalAUMDF$AUM_total_privateandpublic <- finalAUMDF$AUM_private + finalAUMDF$AUM_public

#save	
write.csv(finalAUMDF, file=paste0(forageDir, "AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_final.csv"), row.names=FALSE)



