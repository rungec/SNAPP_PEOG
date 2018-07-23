#This script uses monte-carlo simulation (the Krinsky-Robb method) to estimate marginal effects and confidence intervals of a variable, and then goes on to estimate the impacts on sage grouse habitat.
#We draw n samples of  multivariate normal distribution of the parameter estimates (betas, model coefficients) with means = model estimates and standard deviations drawn from the heteroskedacity-robust variance-covariance matrix. 
#For each draw, we predict the response for each county and a given loss of public grazing land, then backtransform this to get Prc - the probability that range converts to cropland
#The average effect and confidence intervals of the variable on Prc can then be estimated across these n draws
#The marginal effect (for a given variable on Prc) is calculated by approximating the derivative, similar to Stata 
#The response of rangeland rent to loss of public grazing land is a function of the ratio of public to private AUM in a county.
#This script follows on from PerfectEnemyofGood_econometricmodels.r

options(stringsAsFactors=FALSE) # turn off automatic factor coersion
options(scipen=9999) #turn off scientific notation

library(MASS) #for mvnorm
library(parallel) #for mclapply
no_cores=4

#wd <- "C:/Claire/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
wd <- "N:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)

########################
#SET UP DATA
########################
#read data
allCountiesDF <- read.csv("2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static_plus_sgmzs.csv", header=TRUE)
#set categorical variables
	allCountiesDF$Year <- as.factor(allCountiesDF$Year)
	allCountiesDF$Region <- as.factor(allCountiesDF$Region)	
	allCountiesDF$ADMIN_FIPS <- as.factor(allCountiesDF$ADMIN_FIPS)
	names(allCountiesDF)[which(names(allCountiesDF)=="Acres.of.Irrigated.Harvested.Cropland.as.Percent.of.All.Harvested.Cropland.Acreage...2012")] <- "PercentLandIrrigated"
#fill blanks
	allCountiesDF$PercentLandIrrigated[is.na(allCountiesDF$PercentLandIrrigated)] <- 0 #all land
	allCountiesDF$PercentCroplandthatisIrrigated[is.na(allCountiesDF$PercentCroplandthatisIrrigated)] <- 0 #just cropland
#subset out only sage grouse counties & conversion on land capability class 1 to6	
	sgCountiesDF1to6 <- subset(allCountiesDF, SageGrouseCounty==1 & LCC=="LCC1to6")
#drop new years
	sgCountiesDF1to6 <- subset(sgCountiesDF1to6, Year %in% 2008:2012)	
#drop rows with NAs in rent & conversion prop
	sgCountiesDF1to6  <- sgCountiesDF1to6 [complete.cases(sgCountiesDF1to6 [,c("NASS_LandRent_NonIrrigatedCropland","NASS_LandRent_RangelandNEG", "ConversionPropCropRangeV2logZIF0001")]),]
sgCountiesDF1to6 <- droplevels(sgCountiesDF1to6)

#set up prediction data
	sgPredictionDF <- subset(allCountiesDF, SageGrouseCounty==1 & LCC=="LCC1to6")
#fill gaps in range rent with av rent across years
	for(i in which(is.na(sgPredictionDF$NASS_LandRent_RangelandNEG))){
		currADMIN_FIPS <- as.character(sgPredictionDF[i,"ADMIN_FIPS"])
		currsub <- subset(sgPredictionDF, ADMIN_FIPS==currADMIN_FIPS)
		sgPredictionDF$NASS_LandRent_RangelandNEG[i] <- mean(currsub$NASS_LandRent_RangelandNEG, na.rm=TRUE)
	}
#fill gaps in crop rent with av rent across years
	for(i in which(is.na(sgPredictionDF$NASS_LandRent_NonIrrigatedCropland))){
		currADMIN_FIPS <- as.character(sgPredictionDF[i,"ADMIN_FIPS"])
		currsub <- subset(sgPredictionDF, ADMIN_FIPS==currADMIN_FIPS)
		sgPredictionDF$NASS_LandRent_NonIrrigatedCropland[i] <- mean(currsub$NASS_LandRent_NonIrrigatedCropland, na.rm=TRUE)
	}
#fill gaps with surrounding counties
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==46019,"NASS_LandRent_RangelandNEG"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56011, 46063, 30011),"NASS_LandRent_RangelandNEG"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==56023,"NASS_LandRent_RangelandNEG"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56041, 56035, 56039, 16019, 16029, 16007, 49033),"NASS_LandRent_RangelandNEG"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==56037,"NASS_LandRent_RangelandNEG"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56041, 56007, 56035, 56013, 49043, 49009, 08081),"NASS_LandRent_RangelandNEG"])
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==56037,"NASS_LandRent_NonIrrigatedCropland"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56041, 56007, 56023, 56035, 56013, 49043, 49009, 08081),"NASS_LandRent_NonIrrigatedCropland"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==16065,"NASS_LandRent_NonIrrigatedCropland"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(16043, 16051, 16019),"NASS_LandRent_RangelandNEG"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==41025,"NASS_LandRent_NonIrrigatedCropland"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(41037, 41017, 41013, 41023, 41045, 32013,32031),"NASS_LandRent_RangelandNEG"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==56009,"NASS_LandRent_NonIrrigatedCropland"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56025, 56005, 56019, 56045, 56027, 56031, 56001, 56007),"NASS_LandRent_RangelandNEG"])	
sgPredictionDF[sgPredictionDF$ADMIN_FIPS==56017,"NASS_LandRent_NonIrrigatedCropland"] <- mean(sgPredictionDF[sgPredictionDF$ADMIN_FIPS %in% c(56013, 56029, 56043, 56025),"NASS_LandRent_RangelandNEG"])	
	
#drop California - no AUM data for this state
sgPredictionDF <- sgPredictionDF[sgPredictionDF$State!="CA", ]
	
#drop remaining NAs	
sgPredictionDF <- sgPredictionDF[complete.cases(sgPredictionDF[,c("NASS_LandRent_NonIrrigatedCropland","NASS_LandRent_RangelandNEG", "ConversionPropCropRangeV2logZIF0001", "AUM_public", "AUM_private")]),]
sgPredictionDF <- droplevels(sgPredictionDF)

###########################
## DEFINE FUNCTION: HETEROSKEDASCITIY-ROBUST standard error calculation.
###########################
#https://thetarzan.wordpress.com/2011/05/28/heteroskedasticity-robust-and-clustered-standard-errors-in-r/
summaryhrse <- function(model) {
	s <- summary(model)
	X <- model.matrix(model)
	u2 <- residuals(model)^2
	XDX <- 0
	 
	## Here one needs to calculate X'DX. But due to the fact that
	## D is huge (NxN), it is better to do it with a cycle.
	for(i in 1:nrow(X)) {
	XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
	}
	XX1 <- solve(t(X)%*%X) # inverse(X'X)
	varcovar <- XX1 %*% XDX %*% XX1 # Variance calculation (Bread x meat x Bread)
	dfc <- sqrt(nrow(X))/sqrt(nrow(X)-ncol(X))  # degrees of freedom adjustment
	stdh <- dfc*sqrt(diag(varcovar)) # Standard errors of the coefficient estimates are the
	# square roots of the diagonal elements
	 
	t <- model$coefficients/stdh
	p <- 2*pnorm(-abs(t))
	coefresults <- cbind(model$coefficients, stdh, t, p)
	dimnames(coefresults) <- dimnames(s$coefficients)
	results <- list(coeffs=data.frame(coefresults), vcovar=varcovar)
	results
}
###########################
#DEFINE FUNCTION: KRINSKY-ROBB CONFIDENCE INTERVALS
###########################
#A form of montecarlo estimation
#set up the setstep function
setstep <- function(x) {
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}
#set up the krinsky robb function
KrinRobb <- function(chosenmodel, nsim, newdata, variable) {
	vcovmatrix <- summaryhrse(chosenmodel)$vcovar
	estimates <- summaryhrse(chosenmodel)$coeffs[,"Estimate"]

	#generate the  multivariate normal distributed coefficient estimates 
	betas <- mvrnorm(n=nsim, mu=estimates, Sigma=vcovmatrix)
	
	#for each draw of the parameter estimates, estimate marginal effect of variable on Prc
	Prc_estimates <- matrix(unlist(lapply(seq_along(1:nsim), function(i) {
		tmpmodel <- chosenmodel
		tmpmodel$coefficients <- betas[i,]
		
		#set up the data
		d0 <- d1 <- newdata
		d0[[variable]] <- d0[[variable]] - setstep(d0[[variable]])
		d1[[variable]] <- d1[[variable]] + setstep(d1[[variable]])

		P0 <- predict(tmpmodel, newdata = d0)
		P1 <- predict(tmpmodel, newdata = d1)

		#P1Prc <- 1+ZIF - {1 / (1 + exp(P1) )} #conversion prob
		#P0Prc <- 1+ZIF - {1 / (1 + exp(P0) )} #conversion prob

		P1Prc <- (exp(P1)-ZIF)/(1-ZIF+exp(P1)) #prop range to crop
		P0Prc <- (exp(P0)-ZIF)/(1-ZIF+exp(P0)) #prop range to crop		

		#calculate marginal effect of variable on Prc 
		out <- (P1Prc - P0Prc) / (d1[[variable]] - d0[[variable]])
		return(out)
	})), byrow=TRUE, nrow=nsim) #each column represents a county, row i is the marginal effect for i in nsim runs.
	#each row of Prc_stats corresponds to a county
	Prc_stats <- data.frame(marginal_effect=apply(Prc_estimates, 2, FUN=function(x) mean(x)), lowerCI=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[1])), upperCI=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[2])))
	return(Prc_stats)
}

###########################
#SET UP FUNCTION TO ESTIMATE MARGINAL EFFECT FOR 0-100% PUBLIC LAND LOSS
###########################
#calculate the marginal effect of variable for a given loss of public grazing land
marginalFun <- function(chosenmodel, dataset, variable, nsim, predictionyear) {

	marginalEffectsList <- mclapply(0:100, function(percentBLMloss) {
		#set up prediction data based on prediction year control variables
			predictionData <- subset(dataset, Year==predictionyear)
		#convert % BLM loss to proportion BLM retained
			propBLMretained <- (100-percentBLMloss)/100
		#Adjust rangeland rent according to formula newrent=oldrent - oldrent*w*percentBLMloss/100
		#where w is the proportion of private AUM in a county
			predictionData$NASS_LandRent_RangelandNEG <- with(predictionData, NASS_LandRent_RangelandNEG * (AUM_public*propBLMretained + AUM_private)/(AUM_private + AUM_public) )
		#calculate the marginal effect and confidence interval
			ME <- KrinRobb(chosenmodel=chosenmodel, nsim=nsim, newdata=predictionData, variable=variable)
		
	return(data.frame(ADMIN_FIPS=predictionData$ADMIN_FIPS, percentBLMloss=rep(percentBLMloss, nrow(predictionData)), ME))	
	}, mc.cores=no_cores)
		
	marginalEffectsDF <- do.call(rbind, marginalEffectsList)
	return(marginalEffectsDF)
}

###########################
#CALCULATE MARGINAL EFFECT OF VARIABLES BY COUNTY
###########################
#set parameters	
	eps = 1e-7 #the small step across which the instantaneous marginal effect is calcuated
	ZIF=0.0001 #zero inflation factor
	set.seed(1234) #set the seed to get the same random sample each time
	predictionyear=2012

#chosen model
chosenmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated) )

#marginal effect of range rent for 0-100% loss of public grazing, using 2012 control variables
rangeRentMargins <- marginalFun(chosenmod, sgPredictionDF, "NASS_LandRent_RangelandNEG", 10000, 2012)
write.csv(rangeRentMargins, "3 Model output/marginal effects/Marginal_effect_rangelandRent_2012_10000sims.csv", row.names=FALSE)	

#marginal effect of percent cropland in county for 0-100% loss of public grazing, using 2014 control variables
percCropMargins <- marginalFun(chosenmod, sgPredictionDF, "PercentArea_CropStatic", 10000, 2012)
write.csv(percCropMargins, "3 Model output/marginal effects/Marginal_effect_PercentCropland_2012_10000sims.csv", row.names=FALSE)	

#marginal effect of population for 0-100% loss of public grazing, using 2014 control variables
popnMargins <- marginalFun(chosenmod, sgPredictionDF, "Popn", 10000, 2012)		
write.csv(popnMargins, "3 Model output/marginal effects/Marginal_effect_Popn_2012_10000sims.csv", row.names=FALSE)
#####################################		



###########################
#DEFINE FUNCTION: ESTIMATES and CONFIDENCE INTERVALS for change in area converted with change in BLM
###########################
#A form of montecarlo estimation
#in this case, we use a different setstep

KrinRobbArea <- function(chosenmodel, nsim, newdata, percentBLMloss, difference) {
	vcovmatrix <- summaryhrse(chosenmodel)$vcovar
	estimates <- summaryhrse(chosenmodel)$coeffs[,"Estimate"]

	#generate the  multivariate normal distributed coefficient estimates 
	betas <- mvrnorm(n=nsim, mu=estimates, Sigma=vcovmatrix)
	
	#for each draw of the parameter estimates, estimate marginal effect of variable on Prc
	Prc_estimates <- matrix(unlist(lapply(seq_along(1:nsim), function(i) {
		tmpmodel <- chosenmodel
		tmpmodel$coefficients <- betas[i,]
		
		#convert % BLM loss to proportion BLM retained
		propBLMretained <- (100-percentBLMloss)/100
		#set up the data
		#Adjust rangeland rent according to formula newrent=oldrent - oldrent*w*percentBLMloss/100
		#where w is the proportion of private AUM in a county
		d0 <- d1 <- newdata
		d0$NASS_LandRent_RangelandNEG <- newdata$NASS_LandRent_RangelandNEG
		d1$NASS_LandRent_RangelandNEG <- with(newdata, NASS_LandRent_RangelandNEG * (AUM_public*propBLMretained + AUM_private)/(AUM_private + AUM_public) )
		
		P0 <- predict(tmpmodel, newdata = d0)
		P1 <- predict(tmpmodel, newdata = d1)

		#P1Prc <- 1+ZIF - {1 / (1 + exp(P1) )} #conversion prob
		#P0Prc <- 1+ZIF - {1 / (1 + exp(P0) )} #conversion prob

		P1Prc <- (exp(P1)-ZIF)/(1-ZIF+exp(P1)) #prop range to crop
		P0Prc <- (exp(P0)-ZIF)/(1-ZIF+exp(P0)) #prop range to crop		

		#calculate change in conversion prop
		if (difference=="diff") {
			out <- (P1Prc - P0Prc)
		} else if (difference=="percent") {
			out <- 100*(P1Prc - P0Prc)/P0Prc
			out[is.na(out)] <- 0
		} else { 
			out <- c(P1Prc)
		}
		out[out < 0] <- 0
		return(out)
	})), byrow=TRUE, nrow=nsim) #each column represents a county, row i is the marginal effect for i in nsim runs.
	
	if (difference=="percent") {
		Prc_stats <- data.frame(percent_increase_on_background=apply(Prc_estimates, 2, FUN=function(x) mean(x)), percent_increase_on_background_lower=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[1])), percent_increase_on_background_upper=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[2])))
	
	} else {
	#each row of Prc_stats corresponds to a county
	Prc_annualarea_estimates <- apply(Prc_estimates, 1, FUN=function(x) x*predictionData$Area_Ranget0) #annual area lost at t0
	Prc_2050area_estimates <- apply(Prc_estimates, 1, FUN=function(x) {1-(1-x)^(2050-2012)}*predictionData$Area_Ranget0) #area lost 
	
	#summarise 10000 estimates
	Prc_stats <- data.frame(delta_conversion_prop=apply(Prc_estimates, 2, FUN=function(x) mean(x)), lowerCI=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[1])), upperCI=apply(Prc_estimates, 2, FUN=function(x) quantile(x, c(0.25, 0.975)[2])), Area_Ranget0=predictionData$Area_Ranget0, Annual_conversion_ha=apply(Prc_annualarea_estimates, 1, FUN=function(x) mean(x)), Annual_conversion_ha_lower=apply(Prc_annualarea_estimates, 1, FUN=function(x) quantile(x, c(0.25, 0.975)[1])), Annual_conversion_ha_upper=apply(Prc_annualarea_estimates, 1, FUN=function(x) quantile(x, c(0.25, 0.975)[2])),  Converted_by2050_ha=apply(Prc_2050area_estimates, 1, FUN=function(x) mean(x)), Converted_by2050_ha_lower=apply(Prc_2050area_estimates, 1, FUN=function(x) quantile(x, c(0.25, 0.975)[1])), Converted_by2050_ha_upper=apply(Prc_2050area_estimates, 1, FUN=function(x) quantile(x, c(0.25, 0.975)[2])))
	}
	return(Prc_stats)
}

#calculate the additional area converted for a given loss of public grazing land
AreaFun <- function(chosenmodel, dataset, variable, nsim, predictionyear, dif) {

	AreaList <- mclapply(0:100, function(percentBLMloss) {
		#set up prediction data based on prediction year control variables
		predictionData <- subset(dataset, Year==predictionyear)
		#calculate the marginal effect and confidence interval
		ME <- KrinRobbArea(chosenmodel=chosenmodel, nsim=nsim, newdata=predictionData, percentBLMloss=percentBLMloss, difference=dif)
		
	return(data.frame(ADMIN_FIPS=predictionData$ADMIN_FIPS, percentBLMloss=rep(percentBLMloss, nrow(predictionData)), ME))	
	}, mc.cores=no_cores)
	
	AreaDF <- do.call(rbind, AreaList)
	return(AreaDF)
}

###########################
#CALCULATE BACKGROUND CONVERSION
###########################
###Setup parameters
ZIF=0.0001 #zero inflation factor
set.seed(1234) #set the seed to get the same random sample each time
chosenmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated) )

predictionData <- subset(sgPredictionDF, Year==2012)
MEb <- KrinRobbArea(chosenmodel=chosenmod, nsim=10000, newdata=predictionData, percentBLMloss=0, difference="abs")
backgroundconversion <- data.frame(ADMIN_FIPS=predictionData$ADMIN_FIPS, percentBLMloss=rep(0, nrow(predictionData)), MEb)
write.csv(backgroundconversion, "3 Model output/marginal effects/BackgroundConversionRates_byCounty_2012_10000sims.csv", row.names=FALSE)

###########################
#CALCULATE ADDITIONAL AREA CONVERTED WITH CHANGE IN BLM
###########################
###Setup parameters
ZIF=0.0001 #zero inflation factor
set.seed(1234) #set the seed to get the same random sample each time
chosenmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated) )

#marginal effect of range rent on crop area for 0-100% loss of public grazing, using 2012 control variables
rangeRentvAreaDif <- AreaFun(chosenmod, sgPredictionDF, "NASS_LandRent_RangelandNEG", 10000, 2012, dif="diff")
write.csv(rangeRentvAreaDif, "3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_minusbackground.csv", row.names=FALSE)

#marginal effect of range rent on crop area for 0-100% loss of public grazing, using 2012 control variables
rangeRentvAreaDifPerc <- AreaFun(chosenmod, sgPredictionDF, "NASS_LandRent_RangelandNEG", 10000, 2012, dif="percent")
write.csv(rangeRentvAreaDifPerc, "3 Model output/marginal effects/PercentIncreaseonConversion_withgrazingloss_2012_10000sims_minusbackground.csv", row.names=FALSE)

#Effect of range rent on total crop area for 0-100% loss of public grazing, using 2012 control variables
rangeRentvAreaAbs <- AreaFun(chosenmod, sgPredictionDF, "NASS_LandRent_RangelandNEG", 10000, 2012, dif="abs")
write.csv(rangeRentvAreaAbs, "3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_includingbackground.csv", row.names=FALSE)

###########################
# #CALCULATE TIME TO TOTAL CONVERSION
# ###########################
# chosenmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated) )

# #Set up function
# #set up and run function
# counties <- as.character(unique(sgPredictionDF$ADMIN_FIPS))
# timetoconversionDF <- c()
# timetoconversionFun <- function(percentRangeLoss, percentBLMloss){
	# for (currCounty in counties){
	# print(paste0("Starting ", currCounty))
	# #setup loop
	# startPreds <- subset(sgPredictionDF, ADMIN_FIPS==currCounty & Year==2012)
	# startPreds <- droplevels(startPreds)
	# startPreds$Area_Ranget1 <- startPreds$Area_Ranget0
	# #convert % BLM loss to proportion BLM retained
	# propBLMretained <- (100-percentBLMloss)/100
	# startPreds$NASS_LandRent_RangelandNEG <- with(startPreds, NASS_LandRent_RangelandNEG * (AUM_public*propBLMretained + AUM_private)/(AUM_private + AUM_public) )
	# completeconversion <- FALSE
	# counter <- 0
	# #run time loop
	# while(completeconversion==FALSE & counter < 1000){
		# currdf <- startPreds
		# P1 <- predict(chosenmod, newdata = currdf)
		# P1Prc <- (exp(P1)-ZIF)/(1-ZIF+exp(P1))
		# currdf$Area_Ranget1 <-  currdf$Area_Ranget1 - P1Prc*currdf$Area_Ranget1
		# currdf$PercentArea_CropStatic <- 1-(currdf$Area_Ranget1/currdf$TotalRangeorCropAreainLCC_ha)
		# #keep running until 95% of range has been converted	
		# completeconversion <- currdf$Area_Ranget1 < ((100-percentRangeLoss)/100)*currdf$Area_Ranget0
		# counter <- sum(counter, 1)
		# }
		# timetoconversionDF <- rbind(timetoconversionDF,c(currCounty, counter))
		# }
# dimnames(timetoconversionDF)[[2]] <- c("ADMIN_FIPS", "Year_to_Conversion")
# write.csv(timetoconversionDF, paste0(sprintf("3 Model output/marginal effects/TimetoConvert_%spercofRange", percentRangeLoss), sprintf("_%spercBLMloss_byCounty_start2012.csv",  percentBLMloss)), row.names=FALSE)
# return(timetoconversionDF)
# }	

# timetoconversionFun(5,10)
# timetoconversionFun(5,25)
# timetoconversionFun(5,50)
# timetoconversionFun(5,0)
# timetoconversionFun(15,0)
# timetoconversionFun(15,10)
# timetoconversionFun(15,25)
# timetoconversionFun(15,50)
# timetoconversionFun(95,0)	


#################
#ESTIMATE MESIC HABITAT LOSS
#################
#Area of sage grouse mesic habitat in private land in county see PerfectEnemyofGood_createMesicInputs.r
mesic_prop <- read.csv("1 Inputs/Sage grouse/Mesic_habitat_Donnelly/MESIC_proportionConverted_byCounty.csv", header=TRUE)
#mesic_current <- read.csv("1 Inputs/Sage grouse/Mesic_habitat_Donnelly/County_mesic.csv", header=TRUE) #area from 1984 to 2011
rangeRentvAreaAbs <- read.csv("3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_includingbackground.csv", header=TRUE)
rangeRentvAreaDif <- read.csv("3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_minusbackground.csv", header=TRUE)

dataList <- list(rangeRentvAreaAbs, rangeRentvAreaDif)
apps <- list("_includingbackground", "_minusbackground")

lapply(seq_along(dataList), function(i){
	currapp <- apps[[i]]
	#merge mesic habitat areas by county, with area of rangeland converted by 2050 by county
	countyMesicLossDF <- merge(dataList[[i]], mesic_prop, by.x="ADMIN_FIPS", by.y="ADMIN_FIPS", all.x=TRUE)

	#calculate how much mesic habitat remains by 2050
		countyMesicLossDF$Mesic_habitat_remaining_2050_ha <- with(countyMesicLossDF, Area_mesic_ha - Converted_by2050_ha*prop_conversion_on_mesic)
		countyMesicLossDF$Mesic_habitat_remaining_2050_ha_lower <- with(countyMesicLossDF, Area_mesic_ha - Converted_by2050_ha_lower*prop_conversion_on_mesic)
		countyMesicLossDF$Mesic_habitat_remaining_2050_ha_upper <- with(countyMesicLossDF, Area_mesic_ha - Converted_by2050_ha_upper*prop_conversion_on_mesic)
	#replace negative numbers with zero
		countyMesicLossDF[,c("Mesic_habitat_remaining_2050_ha", "Mesic_habitat_remaining_2050_ha_lower", "Mesic_habitat_remaining_2050_ha_upper")] [ countyMesicLossDF[,c("Mesic_habitat_remaining_2050_ha", "Mesic_habitat_remaining_2050_ha_lower", "Mesic_habitat_remaining_2050_ha_upper")] <0 ] <- 0

	#calculate how much mesic habitat is lost by 2050
		countyMesicLossDF$Mesic_habitat_loss_2050_ha <- sapply(1:nrow(countyMesicLossDF), function(x)  min(countyMesicLossDF$Converted_by2050_ha[x]*countyMesicLossDF$prop_conversion_on_mesic[x], countyMesicLossDF$Area_mesic_ha[x]))
		countyMesicLossDF$Mesic_habitat_loss_2050_ha_lower <- sapply(1:nrow(countyMesicLossDF), function(x)  min(countyMesicLossDF$Converted_by2050_ha_lower[x]*countyMesicLossDF$prop_conversion_on_mesic[x], countyMesicLossDF$Area_mesic_ha[x]))
		countyMesicLossDF$Mesic_habitat_loss_2050_ha_upper <- sapply(1:nrow(countyMesicLossDF), function(x)  min(countyMesicLossDF$Converted_by2050_ha_upper[x]*countyMesicLossDF$prop_conversion_on_mesic[x], countyMesicLossDF$Area_mesic_ha[x]))

	#Calculate what % remains by 2050
		countyMesicLossDF$Mesic_habitat_remaining_2050_percent <- with(countyMesicLossDF, 100*Mesic_habitat_remaining_2050_ha/Area_mesic_ha)
		countyMesicLossDF$Mesic_habitat_remaining_2050_percent_lower <- with(countyMesicLossDF, 100*Mesic_habitat_remaining_2050_ha_lower/Area_mesic_ha)
		countyMesicLossDF$Mesic_habitat_remaining_2050_percent_upper <- with(countyMesicLossDF, 100*Mesic_habitat_remaining_2050_ha_upper/Area_mesic_ha)
	#Calculate what % lost by 2050
		countyMesicLossDF$Mesic_habitat_lost_2050_percent <- with(countyMesicLossDF, 100*Converted_by2050_ha*prop_conversion_on_mesic/Area_mesic_ha)
		countyMesicLossDF$Mesic_habitat_lost_2050_percent_lower <- with(countyMesicLossDF, 100*Converted_by2050_ha_lower*prop_conversion_on_mesic/Area_mesic_ha)
		countyMesicLossDF$Mesic_habitat_lost_2050_percent_upper <- with(countyMesicLossDF, 100*Converted_by2050_ha_upper*prop_conversion_on_mesic/Area_mesic_ha)
	#These columns have NAs where no mesic area
	#save
	write.csv(countyMesicLossDF, sprintf("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_model%s.csv", currapp), row.names=FALSE)
	
	#Total loss of mesic habitat
		totalMesicLossDF <- aggregate(countyMesicLossDF[c("Area_Ranget0", "Annual_conversion_ha", 
		"Annual_conversion_ha_upper", "Annual_conversion_ha_lower", "Converted_by2050_ha", 
		"Converted_by2050_ha_upper", "Converted_by2050_ha_lower", "Area_mesic_ha", "Mesic_habitat_remaining_2050_ha", 
		"Mesic_habitat_remaining_2050_ha_lower", "Mesic_habitat_remaining_2050_ha_upper", "Mesic_habitat_loss_2050_ha", "Mesic_habitat_loss_2050_ha_lower", "Mesic_habitat_loss_2050_ha_upper")], by=countyMesicLossDF[c("percentBLMloss")], FUN=sum)
		totalMesicLossDF$Mesic_habitat_loss_2050_percent <- 100 - (100*totalMesicLossDF$Mesic_habitat_remaining_2050_ha/totalMesicLossDF$Area_mesic_ha)
		totalMesicLossDF$Mesic_habitat_loss_2050_percent_lower <- 100 - (100*totalMesicLossDF$Mesic_habitat_remaining_2050_ha_lower/totalMesicLossDF$Area_mesic_ha)
		totalMesicLossDF$Mesic_habitat_loss_2050_percent_upper <- 100 - (100*totalMesicLossDF$Mesic_habitat_remaining_2050_ha_upper/totalMesicLossDF$Area_mesic_ha)
	write.csv(totalMesicLossDF, sprintf("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_total_model%s.csv", currapp), row.names=FALSE)		

	#SENSITIVITY ANALYSIS: ESTIMATE MESIC LOSS - Only for counties where majority first crop is not alfalfa
	#Table made in ArcGIS by tabulate area first_crop_class.tif (first breakout crop, from 2008 to 2012 Lark CDL data) overlapped with county shapefile	
	firstcrop <- read.csv("1 Inputs/Cropland and pasture/First_crop_byCounty_2008to2012_aspercent.csv", header=TRUE)
	alfalfacounties <- firstcrop$ADMIN_FIPS[firstcrop$most_common_breakout_crop=="Alfalfa"]
	
	#in alfalfa counties, set conversion to 0 and remaining habitat to equal area of mesic habitat 
	countyMesicLossDFnoAlf <- countyMesicLossDF
	countyMesicLossDFnoAlf$alfalfacounty <- 0
	countyMesicLossDFnoAlf$alfalfacounty[countyMesicLossDFnoAlf$ADMIN_FIPS %in% alfalfacounties] <- 1

	countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$alfalfacounty==1, c("Annual_conversion_ha", "Annual_conversion_ha_upper", "Annual_conversion_ha_lower", "Converted_by2050_ha", "Converted_by2050_ha_upper", "Converted_by2050_ha_lower", "Mesic_habitat_loss_2050_ha", "Mesic_habitat_loss_2050_ha_lower", "Mesic_habitat_loss_2050_ha_upper")] <- 0
	countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$alfalfacounty==1, c("Mesic_habitat_remaining_2050_ha",  "Mesic_habitat_remaining_2050_ha_lower", "Mesic_habitat_remaining_2050_ha_upper")] <- countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$alfalfacounty==1, "Area_mesic_ha"]
	countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$alfalfacounty==1, c("Mesic_habitat_lost_2050_percent", "Mesic_habitat_lost_2050_percent_lower", 
	"Mesic_habitat_lost_2050_percent_upper")] <- 0
	countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$alfalfacounty==1, c("Mesic_habitat_remaining_2050_percent", 
	"Mesic_habitat_remaining_2050_percent_lower", "Mesic_habitat_remaining_2050_percent_upper")] <- 
	100

	write.csv(countyMesicLossDFnoAlf, sprintf("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_model%s_notAlfalfaCounties.csv", currapp), row.names=FALSE)

	#Total loss of mesic habitat
	totalMesicLossDFnoAlf <- aggregate(countyMesicLossDFnoAlf[c("Area_Ranget0", "Annual_conversion_ha", 
	"Annual_conversion_ha_upper", "Annual_conversion_ha_lower", "Converted_by2050_ha", 
	"Converted_by2050_ha_upper", "Converted_by2050_ha_lower", "Area_mesic_ha", "Mesic_habitat_remaining_2050_ha", 
	"Mesic_habitat_remaining_2050_ha_lower", "Mesic_habitat_remaining_2050_ha_upper")], by=countyMesicLossDFnoAlf[c("percentBLMloss")], FUN=sum)

	totalMesicLossDFnoAlf$Mesic_habitat_loss_2050_percent <- 100 - (100*totalMesicLossDFnoAlf$Mesic_habitat_remaining_2050_ha/totalMesicLossDFnoAlf$Area_mesic_ha) #as percent of mesic habitat in the whole study region, not just the mesic habitat in the non-alfalfa counties
	totalMesicLossDFnoAlf$Mesic_habitat_loss_2050_percent_lower <- 100 - (100*totalMesicLossDFnoAlf$Mesic_habitat_remaining_2050_ha_lower/totalMesicLossDFnoAlf$Area_mesic_ha)
	totalMesicLossDFnoAlf$Mesic_habitat_loss_2050_percent_upper <- 100 - (100*totalMesicLossDFnoAlf$Mesic_habitat_remaining_2050_ha_upper/totalMesicLossDFnoAlf$Area_mesic_ha)

	write.csv(totalMesicLossDFnoAlf, sprintf("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_total_model%s_notAlfalfaCounties.csv", currapp), row.names=FALSE)		
		
})


#################
#SENSITIVITY ANALYSIS: ESTIMATE HABITAT LOSS - extrapolate current conversion rates
#################
library(plyr)
yrid <- list("0812", "1315")
yearlist <- list(2008:2012, 2013:2015)
#Actual conversion, averaged 2008 to 2012
allCountiesDF <- read.csv("2 Data by county/Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables_static_plus_sgmzs.csv", header=TRUE)
#subset out only sage grouse counties & conversion on land capability class 1 to6 and drop CA	
sgCountiesDF1to6 <- subset(allCountiesDF, SageGrouseCounty==1 & LCC=="LCC1to6") 
ZIF=0.0001	
	
lapply(1:2, function(i) {	
	currSub <- subset(sgCountiesDF1to6, Year %in% yearlist[[i]])
	#average conversion across 2008 to 2012
	currSub$Prc <- with(currSub, (exp(ConversionPropCropRangeV2logZIF0001)-ZIF)/(1-ZIF+exp(ConversionPropCropRangeV2logZIF0001)))
	currSub$Actual_annual_converted_ha <- with(currSub, (Prc*Area_Ranget0))
	currSub$Actual_converted_by2050_ha <- with(currSub, {1-(1-Prc)^(2050-2012)}*Area_Ranget0)
	actualConversion <- ddply(currSub, .(ADMIN_FIPS, Region, State), summarize,
					ConversionPropCropRangeV2 = mean(ConversionPropCropRangeV2, na.rm=TRUE),
					Prc = mean(Prc, na.rm=TRUE),
					Area_Ranget0 = round(mean(Area_Ranget0, na.rm=TRUE), 2),
					Actual_annual_converted_ha = round(mean(Actual_annual_converted_ha, na.rm=TRUE), 2),
					Actual_converted_by2050_ha = round(mean(Actual_converted_by2050_ha, na.rm=TRUE), 2)
					)

	#Area of sage grouse mesic habitat in county
	mesic_prop <- read.csv("1 Inputs/Sage grouse/Mesic_habitat_Donnelly/MESIC_proportionConverted_byCounty.csv", header=TRUE)
	#merge mesic habitat areas by county, with area of rangeland converted by 2050 by county
	actualMesicLossDF <- merge(actualConversion, mesic_prop, by.x="ADMIN_FIPS", by.y="ADMIN_FIPS", all.x=TRUE)
	#calculate how much mesic habitat lost each year
	actualMesicLossDF$Actual_annual_mesic_habitat_loss_ha <- with(actualMesicLossDF, Area_mesic_ha*Prc*prop_conversion_on_mesic)
	#calculate how much mesic habitat remains by 2050
	actualMesicLossDF$Actual_mesic_habitat_remaining_2050_ha <- with(actualMesicLossDF, Area_mesic_ha - Actual_converted_by2050_ha*prop_conversion_on_mesic)
	#replace negative numbers with zero
	actualMesicLossDF$Actual_mesic_habitat_remaining_2050_ha[ actualMesicLossDF$Actual_mesic_habitat_remaining_2050_ha <0 ] <- 0
	actualMesicLossDF$Actual_mesic_habitat_remaining_2050_percent <- with(actualMesicLossDF, 100*Actual_mesic_habitat_remaining_2050_ha/Area_mesic_ha)
	actualMesicLossDF$Actual_mesic_habitat_remaining_2050_ha <- with(actualMesicLossDF, 100 - (100*Actual_mesic_habitat_remaining_2050_ha/Area_mesic_ha))
	#save
	write.csv(actualMesicLossDF, sprintf("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_extrapolated%s.csv", yrid[[i]]), row.names=FALSE)
	})


