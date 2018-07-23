#Uses data collated & summarised in PerfectEnemyofGood_mergeData.r

library(MASS) #for boxcox, stepAIC and rlm
#library(glmulti) #for glmulti
library(MuMIn) #for AICc
library(plyr) #for ddply
library(ggplot2)
#library(gridExtra)#for grid.arrange
#library(lme4) #for lmer

#options(stringsAsFactors=FALSE) # turn off automatic factor coersion

#wd <- "Y:/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)

###EDITME

currdate <- "170612_static_MZregion_2008_2012"


########################
###MODEL LARK CDL CONVERSION RATES AGAINST NASS LAND RENT DATA
########################
#NASS land rent data is in $/acre (trialed models with single variable representing deltaRentCropRange, but this is not strictly correct according to the theoretical econometric model, so I ran other models with separate variables for crop & range rents)
##ln(R->C/R->R) ~ C + D(i in State)a(s) + b(t) + d(c)r(c,it) + d(r)r(r,it) + E(it)
##ln conversion prop ~ intercept + fixed effect of state + fixed effect of time + rent crop - rent rangeland + error

#Proportion of rangeland converted (conversionprop) was calculated from data in hectares, but is unitless
#Models using revised conversion proportion #ConversionPropCropRangeV2 (rangeland t0 converted to cropland t1/rangeland t0 that stays rangeland t1)
#ConversionPropCropRangeV2logZIF = log((R->C/R->R)+1)

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

########################
#CHECK CORRELATION between rangeland rent & cropland rent
########################
# with(sgCountiesDF1to6, cor(NASS_LandRent_RangelandNEG, NASS_LandRent_NonIrrigatedCropland))
# #[1] -0.7010497
	# png(filename=paste0(wd, "3 Model output/econometric model/Plot_rangelandrent_vs_croplandrent.png"), width=670, height=670)
		# with(sgCountiesDF1to6, plot(NASS_LandRent_RangelandNEG, NASS_LandRent_NonIrrigatedCropland))
	# dev.off()

########################
#CHECK CORRELATION between rangeland rent & response variable
########################	

	# png(filename=paste0(wd, "3 Model output/econometric model/Plot_rangelandrent_vs_conversion.png"), width=670, height=670)
		# ggplot(sgCountiesDF1to6, aes(x=NASS_LandRent_Rangeland, y=ConversionPropCropRangeV2))+
		# geom_point(size=3)+
		# labs(y="Proportion rangeland converted", x="Rangeland Rent") +
		# ylim(0,0.05)+
		# theme_classic(27)+
		# guides(color=guide_legend(title=NULL))+
		# theme(axis.line.x = element_line(color="grey70"), axis.line.y = element_line(color="grey70"))
	# dev.off()
	# png(filename=paste0(wd, "3 Model output/econometric model/Plot_rangelandrent_vs_conversionZIF0001.png"), width=670, height=670)
		# ggplot(sgCountiesDF1to6, aes(x=NASS_LandRent_Rangeland, y=ConversionPropCropRangeV2logZIF0001))+
		# geom_point(size=3)+
		# labs(y="Proportion rangeland converted", x="Rangeland Rent") +
		# theme_classic(27)+
		# guides(color=guide_legend(title=NULL))+
		# theme(axis.line.x = element_line(color="grey70"), axis.line.y = element_line(color="grey70"))
	# dev.off()


########################
#FULL MODEL, all variables
########################
fullmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Region) + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated + TotalRangeorCropAreainLCC_ha + road_density + Popn + Popn_Chg + PopnChg_Perc + AREA_URBAN_ha + AREAPCT_URBAN + I(PercentArea_CropStatic^2) ))
AICc(fullmod) #1397.218 #839.3781 #1440.207

#########################
#PARSIMONOUS Models
#########################
basemod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + factor(Year) +factor(Region) + NASS_LandRent_NonIrrigatedCropland  + PercentArea_CropStatic + Popn + PercentLandIrrigated + AREA_URBAN_ha + road_density ) )
AICc(basemod) #1416.503 #845.2395 #1445.75
minimalmod <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG  + NASS_LandRent_NonIrrigatedCropland + factor(Year)  + PercentLandIrrigated) )
AICc(minimalmod) #1492.911 #924.0333 #1530.276

########################### ConversionPropCropRangeV2logZIF0001
#SET OF PLAUSIBLE CANDIDATE MODELS
###########################
#Note that are I haven't removed NAs from the dataset sgCountiesDF1to6, all these models will have differing numbers of datapoints, depending on how many NAs are in the variables used in the model.
	#control for area of cropland - most parsimonious model
	mod1 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic ))
	#control for area crop, & percent irrigated
	mod2 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated ))
	#control for area crop, percent irrigated and area of urban
	mod3 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated + AREA_URBAN_ha) )
	#control for area crop, percent irrigated and percent urban
	mod4 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated +  AREAPCT_URBAN) )
	#control for area crop, percent irrigated and road density
	mod5 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated + road_density ) )
	#without state,  control for area crop, percent irrigated, urban, road density and popn
	mod6 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated + Popn + AREA_URBAN_ha) )
	#control for area crop, percent irrigated, area urban and popn change
	mod7 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + PercentArea_CropStatic + PercentLandIrrigated + Popn_Chg + AREA_URBAN_ha) )
	#control for area crop, percent irrigated, state, percent urban and popn change
	mod8 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn_Chg + AREAPCT_URBAN) )
	#control for area crop, percent irrigated, state, area urban
	mod9 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + AREA_URBAN_ha) )
	#full model
	mod10 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha) )
	#full model plus interaction term
	mod10a <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + Popn*AREA_URBAN_ha) )
	#full model plus interaction term
	mod10b <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated) )
	#full model plus saturation fn on crop area
	mod10c <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + I(PercentArea_CropStatic^2)) )
	#without popn
	mod11 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + AREA_URBAN_ha) )
	#without area urban
	mod12 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn ) )
	#without percent irrigated
	mod13 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic  + road_density + Popn + AREA_URBAN_ha) )
	#without road density
	mod14 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentArea_CropStatic + PercentLandIrrigated + Popn + AREA_URBAN_ha) )
	#without percent area crop
	mod15 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) + factor(Region) + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha ) )
	#without year or state
	mod16<- lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland  + PercentArea_CropStatic + Popn + AREA_URBAN_ha + road_density, data=sgCountiesDF1to6)
	mod17<- lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland  + PercentArea_CropStatic + AREA_URBAN_ha, data=sgCountiesDF1to6)
	mod18<- lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland  + PercentArea_CropStatic, data=sgCountiesDF1to6)
	mod19<- lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland  + PercentArea_CropStatic + Popn + AREA_URBAN_ha, data=sgCountiesDF1to6)
	mod20<- lm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + PercentArea_CropStatic + Popn + AREA_URBAN_ha, data=sgCountiesDF1to6)
	mod21<- lm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + PercentArea_CropStatic + Popn + AREA_URBAN_ha + road_density, data=sgCountiesDF1to6)
	mod22<- lm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Year + PercentArea_CropStatic + Popn + AREA_URBAN_ha + road_density, data=sgCountiesDF1to6)
	mod23 <- lm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + PercentArea_CropStatic + Popn + AREA_URBAN_ha + road_density + PercentLandIrrigated, data=sgCountiesDF1to6)
	#full model plus interaction term, no year
	mod24 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Region) +  PercentArea_CropStatic  + road_density + Popn + AREA_URBAN_ha + PercentLandIrrigated+ PercentArea_CropStatic*PercentLandIrrigated) )
#full model plus interaction term, no state
	mod25 <- with(sgCountiesDF1to6, lm(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + factor(Year) +  PercentArea_CropStatic  + road_density + Popn + AREA_URBAN_ha + PercentLandIrrigated+ PercentArea_CropStatic*PercentLandIrrigated) )
	
models <- list(mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4, mod5=mod5, mod6=mod6, mod7=mod7, mod8=mod8, mod9=mod9, mod10=mod10, mod10a=mod10a, mod10b=mod10b,mod24=mod24, mod25=mod25, mod10c=mod10c, mod11=mod11, mod12=mod12, mod13=mod13, mod14=mod14, mod15=mod15, mod16=mod16, mod17=mod17, mod18=mod18, mod19=mod19, mod20=mod20, mod21=mod21, mod22=mod22, mod23=mod23)

# #How well do the models perform
modstats <- cbind(AIC=unlist(lapply(models, function(x) AICc(x))),
BIC=unlist(lapply(models, function(x) BIC(x))))
write.csv(modstats, sprintf("3 Model output/econometric model/AIC_BIC_ofcandidatemodels_%s.csv", currdate))
#save model summaries to a text file
sink(sprintf("3 Model output/econometric model/Summaryofcandidatemodels_ZIF0001_%s.txt", currdate))
print(modstats)
for (i in seq_along(models)){
	print(names(models)[[i]])
	print(summary(models[[i]]))
	print(paste("AIC", names(models)[[i]], round(AICc(models[[i]]), 2), sep=" "))
	print(paste("BIC", names(models)[[i]], round(BIC(models[[i]]), 2), sep=" "))	
	}
sink()

###########################
## HETEROSKEDASCITIY-ROBUST standard error calculation.
###########################
#http://www.statmethods.net/stats/regression.html
# https://www.r-bloggers.com/model-selection-and-multi-model-inference/
############################
#Compare the robust errors with ordinary least squares - if different use robust
#All observations not shown above have a weight of 1. In OLS regression, all cases have a weight of 1. Hence, the more cases in the robust regression that have a weight close to one, the closer the results of the OLS and robust regressions.
# http://statistics.ats.ucla.edu/stat/r/dae/rreg.htm
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
	results <- list(coeffs=coefresults, vcov=varcovar)
	results
}

###########################
#SENSITIVITY ANALYSIS OF COEF for RANGELAND RENT 
###########################
# averaged across candidate models
# #In our case, it makes sense to equally weight the candidate models
#exclude model 10a, as it is computationally singular 
	models <- models[which(names(models)!="mod10a")]
#calculate the heteroskedasticity-robust coefficents and standard errors
	hrsemodels <- lapply(seq_along(models), function(i){
			currestimate <- data.frame(coef=row.names(summaryhrse(models[[i]])$coeffs), summaryhrse(models[[i]])$coeffs, model=rep(names(models)[[i]], nrow(summaryhrse(models[[i]])$coeffs)))
			return(currestimate)
	})
#merge the summaries into a single dataset
	mergedhrseDF <- Reduce(function(...) merge(..., all=TRUE), hrsemodels)
#save
	write.csv(mergedhrseDF, file=sprintf("3 Model output/econometric model/Heteroskedasticityrobust_coefs_forallcandidatemodels_ZIF0001_%s.csv", currdate), row.names=FALSE)

#extract the estimates for the coefficient of rangeland rent
	rangelandrentDF <- mergedhrseDF[mergedhrseDF$coef=="NASS_LandRent_RangelandNEG", ]
	croplandrentDF <- mergedhrseDF[mergedhrseDF$coef=="NASS_LandRent_NonIrrigatedCropland", ]
	rangecoef = mean(rangelandrentDF$Estimate)
	rangecoefmin = min(rangelandrentDF$Estimate)
	rangecoefmax = max(rangelandrentDF$Estimate)
	rangesterror = mean(rangelandrentDF$Std..Error)

#plot the coeficients for rangeland and cropland rent
plotfun <- function(name, x, y, sdev, labelcol){
	plot(x,y, 
		ylim=range(c(y-sdev, y+sdev)),
		main=name,
		xlab="",
		ylab="Estimated coefficent",
		xaxt="n", #surpress x axis
		pch=19,
		cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
	axis(1, at=x, labels=labelcol, las=2)
	# horizontal error bars
	arrows(x, y-sdev, x, y+sdev, length=0.1, angle=90, code=3)
}

png(filename=paste0(wd, sprintf("3 Model output/econometric model/Rent_coefficients_crop_range_modelaverages_ZIF0001_%s.png", currdate)), width=1340, height=670)
par(mfrow=c(1,2))
	#plot rangeland rent coef
	plotfun("Rangeland rent", 1:nrow(rangelandrentDF), rangelandrentDF$Estimate, rangelandrentDF$Std..Error, rangelandrentDF$model)
	#plot cropland rent coef
	plotfun("Cropland rent", 1:nrow(croplandrentDF), croplandrentDF$Estimate, croplandrentDF$Std..Error, croplandrentDF$model)
dev.off()


########################### 
#CHOSEN MODEL - explore different zero inflation transformations
###########################
mod10b<- lm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + Year + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated, data=sgCountiesDF1to6)
mod10ba<- lm(ConversionPropCropRangeV2logZIF1 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + Year + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated, data=sgCountiesDF1to6)
mod10bb<- lm(ConversionPropCropRangeV2logZIF001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + Year + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated, data=sgCountiesDF1to6)
mod10bc<- lm(ConversionPropCropRangeV2logZIFsmall ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + Year + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated, data=sgCountiesDF1to6)

modzifs <- list(z0001=mod10b, z1=mod10ba, z001=mod10bb, zsmall=mod10bc)
#plot models
png(sprintf("3 Model output/econometric model/Summaryofchosenmodel_ZIFexploration_QQplot_%s.png", currdate), height=1240, width=1240)
par(mfrow=c(4,4))
	for (x in seq_along(modzifs)){
		plot(modzifs[[x]], main=names(modzifs[x]) )
}
dev.off()
#write model summary to file
sink(sprintf("3 Model output/econometric model/Summaryofchosenmodel_ZIFexploration_%s.txt", currdate))			
	for (x in seq_along(modzifs)) { 
		print(paste(names(modzifs[x]), "AIC", AICc(modzifs[[x]]), sep=" "))
		print(summary(modzifs[[x]]))
		}
sink()		
#print AICs			
lapply(modzifs, function(x) AICc(x))
lapply(modzifs, function(x) BIC(x))

########################### 
#CHOSEN MODEL - explore hrse
###########################

modrobust <- rlm(ConversionPropCropRangeV2logZIF0001 ~ 1 + NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + Year + PercentArea_CropStatic + PercentLandIrrigated + road_density + Popn + AREA_URBAN_ha + PercentArea_CropStatic*PercentLandIrrigated, data=sgCountiesDF1to6, psi=psi.huber) #for huber robust
hweights <- data.frame(ADMIN_FIPS = sgCountiesDF1to6$ADMIN_FIPS, resid = modrobust$resid, weight = modrobust$w)
	png(filename=paste0(wd, sprintf("3 Model output/econometric model/Heterorobust_datapoint_weights_%s.png", currdate)), width=670, height=670)
		plot(hweights)
	dev.off()
hweights2 <- hweights[order(modrobust$w), ]
hweights2[1:72, ] #71 values are adjusted ie have weight other than 1
summaryhrse(basemod)$coeffs
summaryhrse(minimalmod)$coeffs
summaryhrse(fullmod)$coeffs

sink(sprintf("3 Model output/econometric model/Summaryofchosenmodel_HRSE_%s.txt", currdate))
print(summaryhrse(mod10b)$coeffs)
print(summary(mod10b))
sink()

##########################
#TRIAL MIXED EFFECTS MODELS
##########################	
#test of including control variables as random effects. 
#Benefits: model is not trying to estimate betas (coefficients) for each variable, freeing up degrees of freedom 
#Cons: difficult to justify how to split (cluster) data
library(lme4)

#remove outliers and split data into  2 clusters	
	sgCountiesDF1to6[row.names(sgCountiesDF1to6) %in% c(8808, 27582, 44305, 8616, 16077, 30015, 56011, 32012),1:12]
	nooutliers <- sgCountiesDF1to6[sgCountiesDF1to6$ConversionPropCropRangeV2 < 0.04,]
	ks <- kmeans(nooutliers[,c("ConversionPropCropRangeV2")], 2)
	ksdat <- cbind(nooutliers, ks[1])
	png(filename=paste0(wd, sprintf("3 Model output/econometric model/Plot_mixedeffects_2clustersnooutliers_%s.png", currdate)), width=670, height=670)
	ggplot(ksdat, aes(y=ConversionPropCropRangeV2, x=NASS_LandRent_RangelandNEG, color=cluster))+
	  geom_point()
	dev.off()
	modm5 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat) 
	modm6 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm7 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year +(1|PercentArea_CropStatic), data=ksdat) 
	modm8 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat)
	modm9 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm10 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  Year +(1|PercentArea_CropStatic), data=ksdat)
	sink(sprintf("3 Model output/econometric model/SummaryofMixedEffectsModels_2clustersnooutliers_%s.txt", currdate))
		print(summary(modm5))
		print(summary(modm6))
		print(summary(modm7))
		print(summary(modm8))
		print(summary(modm9))
		print(summary(modm10))
	sink()	

#split data into two clusters, all data points	
	ks <- kmeans(sgCountiesDF1to6[,c("ConversionPropCropRangeV2")], 2)
	ksdat <- cbind(sgCountiesDF1to6, ks[1])
	png(filename=paste0(wd, sprintf("3 Model output/econometric model/Plot_mixedeffects_2clusters_%s.png", currdate)), width=670, height=670)
	ggplot(ksdat, aes(y=ConversionPropCropRangeV2, x=NASS_LandRent_RangelandNEG, color=cluster))+
	  geom_point()
	dev.off()
	modm5 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat) 
	modm6 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm7 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year +(1|PercentArea_CropStatic), data=ksdat) 
	modm8 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat)
	modm9 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm10 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  Year +(1|PercentArea_CropStatic), data=ksdat)
	sink(sprintf("3 Model output/econometric model/SummaryofMixedEffectsModels_2clusters_%s.txt", currdate))
		print(summary(mod5))
		print(summary(mod6))
		print(summary(mod7))
		print(summary(mod8))
		print(summary(mod9))
		print(summary(mod10))
	sink()	

#split data into three clusters, all data points	
	ks <- kmeans(sgCountiesDF1to6[,c("ConversionPropCropRangeV2")], 3)
	ksdat <- cbind(sgCountiesDF1to6, ks[1])
	png(filename=paste0(wd, sprintf("3 Model output/econometric model/Plot_mixedeffects_3clusters_%s.png", currdate)), width=670, height=670)
	ggplot(ksdat, aes(y=ConversionPropCropRangeV2, x=NASS_LandRent_RangelandNEG, color=cluster))+
	  geom_point()
	dev.off()
	modm5 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat) 
	modm6 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm7 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year +(1|PercentArea_CropStatic), data=ksdat) 
	modm8 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + cluster + Year + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=ksdat)
	modm9 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  (1|Year) +(1|PercentArea_CropStatic), data=ksdat) 
	modm10 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  Year +(1|PercentArea_CropStatic), data=ksdat)
	sink(sprintf("3 Model output/econometric model/SummaryofMixedEffectsModels_3clusters_%s.txt", currdate))
		print(summary(modm5))
		print(summary(modm6))
		print(summary(modm7))
		print(summary(modm8))
		print(summary(modm9))
		print(summary(modm10))
	sink()
	
#split data by state, all data points	
	modm5 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=sgCountiesDF1to6) 
	modm6 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region  +(1|PercentArea_CropStatic), data=sgCountiesDF1to6) 
	modm7 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region  + (1|Year) + (1|PercentArea_CropStatic), data=sgCountiesDF1to6) 
	modm8 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region  + (1|Year) + (1|PercentArea_CropStatic) + (1|Popn) + (1|AREAPCT_URBAN), data=sgCountiesDF1to6)
	modm9 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland + Region + (1|PercentArea_CropStatic) + (1|Popn) , data=sgCountiesDF1to6) 
	modm10 <- lmer(ConversionPropCropRangeV2logZIF0001 ~ NASS_LandRent_RangelandNEG + NASS_LandRent_NonIrrigatedCropland +  Region +(1|PercentArea_CropStatic) + (1|AREAPCT_URBAN), data=sgCountiesDF1to6)
	sink(sprintf("3 Model output/econometric model/SummaryofMixedEffectsModels_stateclusters_%s.txt", currdate))
		print(summary(modm5))
		print(summary(modm6))
		print(summary(modm7))
		print(summary(modm8))
		print(summary(modm9))
		print(summary(modm10))
	sink()
#the predictions of these models are very highly dependent on the clusters chosen, and it is difficult to find a justifiable criteria for splitting the data	
#for this reason we continue with a fixed effects model rather than mixed effects	

#AT THIS POINT, SWITCH TO SCRIPT PerfectEnemyofGood_model_marginaleffects.R
#CALCULATIONS BELOW FOR PLOTS FOR MARCH 2017 WORKING GROUP MEETING
###########################
#Relationship between AUM & profit
###########################
#data from Taylor, D.T., Coupal, R.H. & Foulke, T. (2005). The Economic Impact of Federal Grazing on the Economy of Park County, Wyoming. University of Wyoming Department of Agricultural and Applied Economics.
#for Park County, Wyoming. We use wyoming as it is representative of the study region; Table 3 p14
#we ignore the impacts of debt and off-farm income on farm profit

RanchProfit <- c(30347, 27506, 22823, 13674, 2463)
RanchProfitlosspercent <- 100*((30347-RanchProfit)/30347)
AUMpercentloss <- c(0, 10, 25, 50, 100)
#totalAUM <- c(762, 725, 666, 559, 352) #AUM used?
broodCows <- c(504, 480, 441, 369, 229) #number of cattle they can keep on the farm
#AUMSforRanch <- c(State=538, BLM=1882, USFS=1883, Private=500, Deeded=1075)

modelRanchProfitvsAUM <- lm(RanchProfitlosspercent ~ 0 + AUMpercentloss) #intercept through zero

sink("3 Model output/ModelSummary_AUM_vs_profit.txt")
summary(modelRanchProfitvsAUM)
sink()

#plot(modelRanchProfitvsAUM)

modelRanchProfitvsCows <- lm(RanchProfitlosspercent ~ broodCows)
summary(modelRanchProfitvsCows)
#plot(modelRanchProfitvsCows)

#modelRanchProfitvsCows2 <- lm(RanchProfitpercent ~ broodCows + I(broodCows^2))
#summary(modelRanchProfitvsCows2)
#plot(modelRanchProfitvsCows2)
#modelAUMvsCows <- lm(AUMpercent ~ broodCows)

png(filename=paste0(wd, "3 Model output/RanchProfit_vs_AUM_numCows.png"), width=1340, height=670)
par(mfrow=c(1,2))
	plot(RanchProfitlosspercent ~ AUMpercentloss, xlab="Percent of AUM lost", ylab="Ranch profit lost, as percent of maximum", pch=19)
	abline(modelRanchProfitvsAUM)
	plot(RanchProfitlosspercent ~ broodCows, xlab="Number of brood cows retained", ylab="Ranch profit lost, as percent of maximum", pch=19)
	abline(modelRanchProfitvsCows)
dev.off()
ranchDF <- data.frame(RanchProfit, RanchProfitlosspercent, AUMpercentloss, broodCows)
png(filename=paste0(wd, "3 Model output/RanchProfit_vs_AUM_numCows_ppt.png"), width=1340, height=670)
	p1<- ggplot(ranchDF, aes(x=AUMpercentloss, y=RanchProfitlosspercent))+
		geom_point(size=3)+
		geom_abline(intercept=0, slope=coef(modelRanchProfitvsAUM)[[1]])+
		labs(y="Ranch profit lost, % of max", x="Percent of AUM lost") +
		theme_classic(27)+
		guides(color=guide_legend(title=NULL))+
		theme(axis.line.x = element_line(color="grey70"), axis.line.y = element_line(color="grey70"))
	p2 <- ggplot(ranchDF, aes(x=broodCows, y=RanchProfitlosspercent))+
		geom_point(size=3)+
		geom_abline(intercept=coef(modelRanchProfitvsCows)[[1]], slope=coef(modelRanchProfitvsCows)[[2]])+
		labs(y="Ranch profit lost, % of max", x="Number of brood cows retained") +
		theme_classic(27)+
		guides(color=guide_legend(title=NULL))+
		theme(axis.line.x = element_line(color="grey70"), axis.line.y = element_line(color="grey70"))
	grid.arrange(p1,p2, ncol=2)
dev.off()

###########################
#Calculate conversion rates under changes to AUM on public land
###########################
#we use the relationship above to calculate how much profit is lost from each county with loss of AUM
#We make the simplifying assumption that rangeland rent is directly related to AUM
#remove NAs & subset for 2012
sgCountiesDF1to6NoNAs <- sgCountiesDF1to6[complete.cases(sgCountiesDF1to6[,c("ADMIN_FIPS", "Year", "State", "StableArea_noncrop", "TotalRangeorCropAreainLCC_ha", "LCC", "Area_Ranget0_to_Cropt1", "Area_Ranget0_to_Ranget1", "PercentArea_CropStatic", "ConversionPropCropRangeV2logZIF0001", "NASS_LandRent_Rangeland", "NASS_LandRent_NonIrrigatedCropland",  "NASS_LandRent_RangelandNEG", "SageGrouseCounty", "State", "road_density", "Popn", "Popn_Chg", "AREA_URBAN_ha", "AREAPCT_URBAN", "PercentLandIrrigated")]),]
sgCountiesDF1to6NoNAs <- sgCountiesDF1to6NoNAs[sgCountiesDF1to6NoNAs$Year==2012, ] #select only data for 2012
sgCountiesDF1to6NoNAs <- sgCountiesDF1to6[sgCountiesDF1to6$Year==2012, ] #select only data for 2012
sgCountiesDF1to6NoNAs <- droplevels(sgCountiesDF1to6NoNAs)

backgroundConversion <- mean(100*(sgCountiesDF1to6NoNAs$ConversionPropCropRangeV2))

predictFun <- function(AUMloss, modelID, baseData, ZIFvalue){
	
	#new rent, as fn of how much AUM is lost
	currMod <- models[[modelID]]
	predictionData <- baseData
	predictionData$NASS_LandRent_RangelandNEG <- baseData$NASS_LandRent_RangelandNEG - (baseData$NASS_LandRent_RangelandNEG*coef(modelRanchProfitvsAUM)[[1]][1]*AUMloss/100)

	#predict conversion rate for each model, using new rent
	preds <- predict(currMod, predictionData, interval='predict')
	predsCurrent <- predict(currMod, baseData, interval='predict')
	TotalRangeArea <- predictionData$Area_Ranget0_to_Ranget1
	###THIS BACKCALCULATION IS WRONG!!!
	#RangeConvertedtoCrop_ha <- (exp(preds[,"fit"]-ZIFvalue)*predictionData$Area_Ranget0_to_Ranget1)
	#RangeConvertedtoCrop_min <- (exp(preds[,"lwr"]-ZIFvalue)*predictionData$Area_Ranget0_to_Ranget1)
	#RangeConvertedtoCrop_max <- (exp(preds[,"upr"]-ZIFvalue)*predictionData$Area_Ranget0_to_Ranget1)
	NewConversion_percentofrange <- 100*(RangeConvertedtoCrop_ha / TotalRangeArea)
	NewConversion_percentofrange_min <- 100*(RangeConvertedtoCrop_min / TotalRangeArea)
	NewConversion_percentofrange_max <- 100*(RangeConvertedtoCrop_max / TotalRangeArea)
	#current conversion rate (from data)
	CurrentConversion_ha <- baseData$Area_Ranget0_to_Cropt1
	CurrentConversion_percentofrange <- 100*(baseData$ConversionPropCropRangeV2)
	#predicted current conversion rate (from model)
	#CurrentConversion_ha_predicted <- (exp(predsCurrent[,"fit"]-ZIFvalue)*predictionData$Area_Ranget0_to_Ranget1)
	CurrentConversion_percentofrange_predicted <- 100*(CurrentConversion_ha_predicted/ TotalRangeArea)
	#difference with background
	AdditionalAreaConverted_ha_byCounty <- RangeConvertedtoCrop_ha - CurrentConversion_ha_predicted
	AdditionalAreaConverted_ha_byCounty_min <- RangeConvertedtoCrop_min - CurrentConversion_ha_predicted
	AdditionalAreaConverted_ha_byCounty_max <- RangeConvertedtoCrop_max - CurrentConversion_ha_predicted

	newConversionbyCounty <- data.frame(model=rep(names(models)[modelID], nrow(predictionData)), 
									county=baseData$ADMIN_FIPS,
									oldRent=baseData$NASS_LandRent_RangelandNEG,
									newRent=predictionData$NASS_LandRent_RangelandNEG,
									AUMloss=rep(AUMloss, nrow(predictionData)),
									predsCurrent,
									preds,
									TotalRangeArea,
									RangeConvertedtoCrop_ha, 
									RangeConvertedtoCrop_min, 
									RangeConvertedtoCrop_max, 
									CurrentConversion_ha,
									CurrentConversion_ha_predicted,
									CurrentConversion_percentofrange,
									CurrentConversion_percentofrange_predicted,
									NewConversion_percentofrange,
									NewConversion_percentofrange_min,
									NewConversion_percentofrange_max,
									AdditionalAreaConverted_ha_byCounty,
									AdditionalAreaConverted_ha_byCounty_min,
									AdditionalAreaConverted_ha_byCounty_max)

	return(newConversionbyCounty)
}

#Summarise across all counties by model & %AUM loss
totalfun <- function(modelresults){
	modeltotals <- ddply(modelresults, c("model", "AUMloss"), summarise, 
			RangeConvertedtoCrop_total_ha = sum(RangeConvertedtoCrop_ha, na.rm=TRUE),
			RangeConvertedtoCrop_total_ha_min = sum(RangeConvertedtoCrop_min, na.rm=TRUE),
			RangeConvertedtoCrop_total_ha_max = sum(RangeConvertedtoCrop_max, na.rm=TRUE),
			RangeConvertedtoCrop_total_aspercent = 100*sum(RangeConvertedtoCrop_ha, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			RangeConvertedtoCrop_total_aspercent_min = 100*sum(RangeConvertedtoCrop_min, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			RangeConvertedtoCrop_total_aspercent_max = 100*sum(RangeConvertedtoCrop_max, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			NewConversion_mean_acrossCountyYears = mean(NewConversion_percentofrange, na.rm=TRUE),
			CurrentConversion_mean_acrossCountyYears = mean(CurrentConversion_percentofrange, na.rm=TRUE),
			CurrentConversion_mean_predicted_acrossCountyYears = mean(CurrentConversion_percentofrange_predicted, na.rm=TRUE),
			CurrentConversion_total_ha = sum(CurrentConversion_ha, na.rm=TRUE),
			CurrentConversion_total_predicted_ha = sum(CurrentConversion_ha_predicted, na.rm=TRUE),
			CurrentConversion_total_percentofrange = 100*sum(CurrentConversion_ha, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			CurrentConversion_total_predicted_percentofrange = 100*sum(CurrentConversion_ha_predicted, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			AdditionalAreaConverted_ha_total = sum(AdditionalAreaConverted_ha_byCounty, na.rm=TRUE),
			AdditionalAreaConverted_ha_total_min = sum(AdditionalAreaConverted_ha_byCounty_min),
			AdditionalAreaConverted_ha_total_max = sum(AdditionalAreaConverted_ha_byCounty_max),
			AdditionalPercentofRangeConverted = 100*sum(AdditionalAreaConverted_ha_byCounty, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			AdditionalPercentofRangeConverted_min = 100*sum(AdditionalAreaConverted_ha_byCounty_min, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			AdditionalPercentofRangeConverted_max = 100*sum(AdditionalAreaConverted_ha_byCounty_max, na.rm=TRUE)/sum(TotalRangeArea, na.rm=TRUE),
			ChangeinConversion_total_aspercent = 100*mean((NewConversion_percentofrange -  CurrentConversion_percentofrange_predicted)/CurrentConversion_percentofrange_predicted, na.rm=TRUE),
			ChangeinConversion_total_aspercent_sd = 100*sd((NewConversion_percentofrange -  CurrentConversion_percentofrange_predicted)/CurrentConversion_percentofrange_predicted, na.rm=TRUE)			
	)
	return(modeltotals)
	}

#predict conversion rates for all models and loss of aum from 0 to 100%
a <- lapply(seq_along(models), function(x){
		b <- lapply(0:100, predictFun, modelID=x, baseData=sgCountiesDF1to6NoNAs, ZIFvalue=0.0001)
	btotal <- do.call(rbind, b)
	return(btotal)
})
allmods <- do.call(rbind, a)
#oup <- mapply(predictFun, 0:100, seq_along(models))
#save
write.csv(allmods, "3 Model output/Modelpredictions_allModels_AUMloss0to100_byCounty_ZIF0001.csv", row.names=FALSE)

totals <- totalfun(allmods)
#save
write.csv(totals, "3 Model output/Modelpredictions_allModels_AUMloss0to100_totals_ZIF0001.csv", row.names=FALSE)

#predict conversion rates for all models and loss of aum from 0 to 100%
a_ZIF001 <- lapply(seq_along(models_ZIF001), function(x){
		b <- lapply(0:100, predictFun, modelID=x, baseData=sgCountiesDF1to6NoNAs, ZIFvalue=0.001)
	btotal <- do.call(rbind, b)
	return(btotal)
})
allmods_ZIF001 <- do.call(rbind, a_ZIF001)
#save
write.csv(allmods_ZIF001, "3 Model output/Modelpredictions_allModels_AUMloss0to100_byCounty_ZIF001.csv", row.names=FALSE)

totals_ZIF001 <- totalfun(allmods_ZIF001)
#save
write.csv(totals_ZIF001, "3 Model output/Modelpredictions_allModels_AUMloss0to100_totals_ZIF001.csv", row.names=FALSE)

sgCountiesDF1to6NoNAsnooutliers <- subset(sgCountiesDF1to6NoNAs, !ADMIN_FIPS %in% c(16077, 30015, 56011, 31021))
sgCountiesDF1to6NoNAsnooutliers <- droplevels(sgCountiesDF1to6NoNAsnooutliers)
a_nooutliers <- lapply(seq_along(models_nooutliers), function(x){
		b <- lapply(0:100, predictFun, modelID=x, baseData=sgCountiesDF1to6NoNAsnooutliers, ZIFvalue=0.0001)
	btotal <- do.call(rbind, b)
	return(btotal)
})
allmods_nooutliers <- do.call(rbind, a_nooutliers)
#save
write.csv(allmods_nooutliers, "3 Model output/Modelpredictions_allModels_AUMloss0to100_byCounty_nooutliers.csv", row.names=FALSE)

totals_nooutliers <- totalfun(allmods_nooutliers)
#save
write.csv(totals_nooutliers, "3 Model output/Modelpredictions_allModels_AUMloss0to100_totals_nooutliers.csv", row.names=FALSE)

	
##############################
#PLOT new conversion rates
##############################

png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionPercent_bymodel.png"), width=1340, height=1340)
d <- ggplot(totals, aes(x=AUMloss, y=AdditionalPercentofRangeConverted)) +
	geom_line()+
	geom_ribbon(aes(ymin=AdditionalPercentofRangeConverted_min, ymax=AdditionalPercentofRangeConverted_max), alpha=0.1) +
	#coord_cartesian(ylim=c(0,8))+
	labs(y="Increase on annual rate of rangeland conversion to cropland (% of total rangeland)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
d + facet_wrap(~model)	
dev.off()

png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionArea_bymodel.png"), width=1340, height=1340)
d <- ggplot(totals, aes(x=AUMloss, y=AdditionalAreaConverted_ha_total)) +
	geom_line()+
	geom_ribbon(aes(ymin=AdditionalAreaConverted_ha_total_min, ymax=AdditionalAreaConverted_ha_total_max), alpha=0.1) +
	#coord_cartesian(ylim=c(0,8))+
	labs(y="Increase in rangeland conversion to cropland (ha)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
d + facet_wrap(~model)	
dev.off()

png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionPercent.png"), width=670, height=670)
ggplot(totals, aes(x=AUMloss, y=AdditionalPercentofRangeConverted, group=model)) +
	geom_line()+
	geom_ribbon(aes(ymin=AdditionalPercentofRangeConverted_min, ymax=AdditionalPercentofRangeConverted_max), alpha=0.05) +
	#coord_cartesian(ylim=c(0,8))+
	labs(y="Increase on annual conversion rate (ha new cropland/ha rangeland)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
	dev.off()

#change in conversion rate	
png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionPercent_change.png"), width=670, height=670)
ggplot(totals, aes(x=AUMloss, y=ChangeinConversion_total_aspercent, group=model)) +
	geom_line()+
	geom_ribbon(aes(ymin=(ChangeinConversion_total_aspercent-ChangeinConversion_total_aspercent_sd), ymax=(ChangeinConversion_total_aspercent+ChangeinConversion_total_aspercent_sd)), alpha=0.1) +
	#coord_cartesian(ylim=c(0,8))+
	labs(y="Increase on background conversion of rangeland to cropland (%)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
	dev.off()

#change in conversion rate	by model
png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionPercent_change_bymodel.png"), width=1340, height=1340)
d <- ggplot(totals, aes(x=AUMloss, y=ChangeinConversion_total_aspercent, group=model)) +
	geom_line()+
	geom_ribbon(aes(ymin=(ChangeinConversion_total_aspercent-ChangeinConversion_total_aspercent_sd), ymax=(ChangeinConversion_total_aspercent+ChangeinConversion_total_aspercent_sd)), alpha=0.1) +
	#coord_cartesian(ylim=c(0,8))+
	labs(y="Increase on background conversion of rangeland to cropland (%)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
d+ facet_wrap(~model)
	dev.off()
	
png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionArea.png"), width=670, height=670)
ggplot(totals, aes(x=AUMloss, y=AdditionalAreaConverted_ha_total, group=AUMloss)) +
	geom_boxplot()+
	#geom_ribbon(aes(ymin=ChangeinConversion_total_aspercent_min, ymax=ChangeinConversion_total_aspercent_max), alpha=0.1) +
	#coord_cartesian(ylim=c(0,5))+
	labs(y="Additional area rangeland converted to cropland (ha)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
	dev.off()
	
png(filename=paste0(wd, "3 Model output/ModelPredictions_conversionPercent_boxplot.png"), width=670, height=670)
ggplot(totals, aes(x=AUMloss, y=AdditionalPercentofRangeConverted, group=AUMloss)) +
	geom_boxplot()+
	#geom_ribbon(aes(ymin=ChangeinConversion_total_aspercent_min, ymax=ChangeinConversion_total_aspercent_max), alpha=0.1) +
	#coord_cartesian(ylim=c(0,5))+
	labs(y="Increase on annual conversion rate (ha new cropland/ha rangeland)", x="Loss of AUM (%)") +
	theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
dev.off()
	
#predicted rangetocrop against actual rangetocrop
subdat <- subset(allmods, AUMloss==0)
png(filename=paste0(wd, "3 Model output/ModelPredictions_AreaConverted_actualvspredicted.png"), width=670, height=670)
ggplot(subdat, aes(y=RangeConvertedtoCrop_ha, x=CurrentConversion_ha, group=model))+
		geom_abline(slope=1, intercept=0, linetype=2)+
		geom_point(aes(color=model))+
		labs(x="Actual hectares rangeland converted to cropland", y="Predicted hectares rangeland converted to cropland") +
		theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
dev.off()

png(filename=paste0(wd, "3 Model output/ModelPredictions_ConversionRate_actualvspredicted.png"), width=670, height=670)
ggplot(subdat, aes(y=NewConversion_percentofrange, x=CurrentConversion_percentofrange, group=model))+
		geom_abline(slope=1, intercept=0, linetype=2)+
		geom_point(aes(color=model))+
		labs(x="Actual proportion rangeland converted", y="Predicted proportion rangeland converted") +
		theme_classic(27)+
	guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
dev.off()

##END
######################
