#Make plots for the PEOG paper

library(rgeos)
library(maptools)
library(sf) #for shps
library(ggplot2)
library(tmap) #spatial equivalent of ggplot
library(plyr) #aggregate
library(ggthemes) #for pretty themes
library(RColorBrewer)
library(grid) #for inset map 
library(gridExtra) #for cornerlabels grob
#library(parallel)  
#no_cores <- 4

wd <- "C:/Claire/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
#wd <- "/home/runge/Data/NCEAS_Postdoc/P1 Sage Grouse/Analysis/1_PEOG_grazing_economics/"
setwd(wd)

options(scipen=9999) #turn off scientific notation

#barplot(rep(1,8), yaxt="n", col=palette())
divergingCols <- c("#FDE428", "#FFF6B3", "#A1D0E4", "#667BBE", "#304999")
monoCols1 <- c("#36213E", "#554971", "#63768D", "#8AC6D0", "#B8F3FF")
monoCols2 <- c("#A1D0E4","#628CD6", "#494592", "#66095A", "#9C0C63", "#340027")
monoCols3 <- c("grey95", "#A1D0E4","#628CD6", "#494592", "#66095A", "#9C0C63", "#340027")
monoCols4 <- c("#D3DCD8", "#A1D0E4","#628CD6", "#534DC1", "#221F60")
divCols2<- c("#D3DCD8", "#CDE7BE",  "#D3E298", "#ECDD7B",  "#CE8147", "#561D25"  )
divCols<- c("#D3DCD8", "#6A98AB", "#CDE7BE",  "#ECDD7B",  "#CE8147", "#561D25"  )
seqCols<- c("#D3DCD8", "#95D4B3", "#77C69B", "#6A98AB", "#396B8E", "#647886")
seqCols2 <- c("#D3DCD8", "#FED98E", "#FE9929", "#D95F0E", "#993404")


########################
#SET UP DATA
########################
#load shapefiles
sgCountiesShp <- st_read("1 Inputs/Boundaries/SGCounties/SGCounties_overlapping_PACSandBreeding_NAD83alb.shp")
usStatesShp <- st_read("1 Inputs/Boundaries/USGS_State_boundaries_USContinental_smooth.shp")
publicShp <- st_read("1 Inputs/Land Ownership and Protected Areas/SGCounties/SGCounties_public_smooth_neg1kmbuffer.shp")

#marginal effect of range rent for 0-100% loss of public grazing, using 2012 control variables
rangeRentMargins <- read.csv("3 Model output/marginal effects/Marginal_effect_rangelandRent_2012_10000sims.csv", header=TRUE)	

#marginal effect of percent cropland in county for 0-100% loss of public grazing, using 2014 control variables
percCropMargins <- read.csv("3 Model output/marginal effects/Marginal_effect_PercentCropland_2012_10000sims.csv", header=TRUE)	

#marginal effect of population for 0-100% loss of public grazing, using 2014 control variables
popnMargins <- read.csv("3 Model output/marginal effects/Marginal_effect_Popn_2012_10000sims.csv", header=TRUE)

#impacts on mesic habitat, summed over landscape
totalMesicLossDF <- read.csv("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_total_model_includingbackground.csv", header=TRUE)
#impacts on mesic habitat, summed over counties where alfalfa is not the breakout crop
totalMesicLossDFnoAlf <- read.csv("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_total_model_includingbackground_notAlfalfaCounties.csv", header=TRUE)

#impacts on mesic habitat by county	
countyMesicLossDF <- read.csv("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_model_includingbackground.csv", header=TRUE)#impacts on mesic habitat by county, excluding counties where alfalfa is the breakout crop	
countyMesicLossDF[countyMesicLossDF$Area_mesic_ha==0, c("Mesic_habitat_remaining_2050_ha", "Mesic_habitat_loss_2050_ha")] <- NA
countyMesicLossDFnoAlf <- read.csv("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_model_includingbackground_notAlfalfaCounties.csv", header=TRUE)
countyMesicLossDFnoAlf[countyMesicLossDFnoAlf$Area_mesic_ha==0, c("Mesic_habitat_remaining_2050_ha", "Mesic_habitat_loss_2050_ha")] <- NA

#impacts on mesic habitat by county, excluding background
countyMesicLossDFdif <- read.csv("3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_model_minusbackground.csv", header=TRUE)

#background conversion rates
backgroundconversion <- read.csv("3 Model output/marginal effects/BackgroundConversionRates_byCounty_2012_10000sims.csv", header=TRUE)

#conversion rates for change to BLM policy
#values with background conversion already subtracted
rangeRentvAreaDif <- read.csv("3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_minusbackground.csv", header=TRUE)
#percent values with background conversion already subtracted
rangeRentvAreaDifPerc <- read.csv("3 Model output/marginal effects/PercentIncreaseonConversion_withgrazingloss_2012_10000sims_minusbackground.csv", header=TRUE)

#total conversion, including background conversion
rangeRentvAreaAbs <- read.csv("3 Model output/marginal effects/AdditionalAreaConverted_withgrazingloss_2012_10000sims_includingbackground.csv", header=TRUE)

actualConversion <- "3 Model output/sage grouse impacts/Area_mesic_habitat_lost_byCounty_extrapolated"	

areagrassland <- read.csv("1 Inputs/LULC/Area_grassorsagebrush_bycounty.csv", header=TRUE)			

#################
#SET UP MAP FUNCTION
#################
#map function

mapfun <- function(dataset, currpercent, currfile, fileext, inset, cornerlabel, ...){
	#select data
	if(is.numeric(currpercent)==TRUE){
	currData <- subset(dataset, percentBLMloss==currpercent)
	} else {
	currData <- dataset
	}
	#merge data with .shp
	sgShp <- merge(sgCountiesShp, currData, by.x='fips', by.y='ADMIN_FIPS')
	
	#inset map
	insetlocatormap <- tm_shape(usStatesShp) + 
					tm_polygons(col="grey70", border.col="grey70") +
				tm_shape(sgCountiesShp)+	
					tm_polygons(col="grey30", border.col="grey30")+
					tm_layout(frame=FALSE)

	#plot map
	currmap <- tm_shape(sgShp) +
		tm_fill(..., auto.palette.mapping=FALSE, contrast=1)+
		tm_style_white()+
		tm_layout(title=cornerlabel, title.size=1.6, inner.margins=c(.08,.12, .10, .1), outer.margins=c(0,0,0,0), frame=TRUE,
		#tm_layout(title=cornerlabel, title.size=1.6, inner.margins=c(.1,.12, .02, .1), outer.margins=c(0,0,0,0), frame=TRUE,
				legend.title.size = 1.2,legend.text.size = 0.9, legend.position = c("right","bottom"), legend.bg.color = "white", legend.bg.alpha = 1, asp=1.25) +
		tm_shape(publicShp) + 
			tm_fill(col="white")+
		tm_shape(usStatesShp) +
			tm_borders(lwd=1, col = "grey70", alpha = .5)	+	
			tm_shape(sgShp) +
			tm_borders(lwd=0.5, col = "grey40", alpha=0.5) #county boundaries	

			#tm_text("STATE", col="grey50")
		
	#save_tmap
	outfile <- sprintf(currfile, currpercent)
		if(inset==FALSE) {		  
		save_tmap(currmap, filename=paste0(outfile, fileext), outer.margins=0.01, width=7)
	} else if(inset==TRUE) {
		save_tmap(currmap + tm_layout(inner.margins = c(.1,.17, .02, .12)), filename=paste0(outfile, "_inset", fileext), insets_tm=insetlocatormap, insets_vp=viewport(x=0.61, y=0.04, width=0.17, height=0.17, just=c("left", "bottom")), outer.margins=0.01, width=7)
	}
}	


#################
#TABLE 1: Area of habitat remaining after 20yrs for 0, 10, 25, 50% BLM loss
#################
#Study region totals
#background conversion rate
#new conversion rate
#increase on background
#Area converted by 2050
#Area of sage grouse mesic habitat 
#Area of sage grouse mesic habitat affected
#%sage grouse mesic habitat affected
#averaged across counties
#these numbers exclude counties in CA that don't have AUM data
#Area of grassland, whole region = 1,098,786km2
#Area of grassland, private = 428,933km2

Table1 <- lapply(list(0,10,25,50,100), FUN=function(currperc) {
	rangeRentvAreaAbs <- merge(rangeRentvAreaAbs, areagrassland, by="ADMIN_FIPS", all.x=TRUE)
	rangeRentvAreaDif <- merge(rangeRentvAreaDif, areagrassland, by="ADMIN_FIPS", all.x=TRUE)
	rangeRentvAreaDifPerc <- merge(rangeRentvAreaDifPerc, areagrassland, by="ADMIN_FIPS", all.x=TRUE)
	currRangeSub <- rangeRentvAreaAbs[rangeRentvAreaAbs$percentBLMloss==currperc, ]
	currRangeDif <- rangeRentvAreaDif[rangeRentvAreaDif$percentBLMloss==currperc, ]
	currRangeDifPerc <- rangeRentvAreaDifPerc[rangeRentvAreaDifPerc$percentBLMloss==currperc, ]
	currTotalSub <- totalMesicLossDF[totalMesicLossDF$percentBLMloss==currperc,]
	currTotalSubnoAlf <- totalMesicLossDFnoAlf[totalMesicLossDFnoAlf$percentBLMloss==currperc,]
	currTableA <- with(currRangeSub, data.frame(
		conversion_rate = mean(delta_conversion_prop),
		conversion_min = min(delta_conversion_prop),
		conversion_max = max(delta_conversion_prop),
		conversion_lowerCI = mean(lowerCI),
		conversion_lowerCI_min = min(lowerCI),
		conversion_lowerCI_max = max(lowerCI),
		conversion_upperCI = mean(upperCI),
		conversion_upperCI_min = min(upperCI),
		conversion_upperCI_max = max(upperCI),
		conversion_upperCI_max = max(upperCI),
		annual_converted_ha = sum(Annual_conversion_ha),
		annual_converted_ha_lowerCI = sum(Annual_conversion_ha_lower),
		annual_converted_ha_upperCI = sum(Annual_conversion_ha_upper),	
		convertedby2050 = sum(Converted_by2050_ha),
		convertedby2050_lowerCI = sum(Converted_by2050_ha_lower),
		convertedby2050_upperCI = sum(Converted_by2050_ha_upper),
		area_grassland=sum(AREA_GRASS_HA),
		convertedby2050_percent = 100*sum(Converted_by2050_ha)/sum(AREA_GRASS_HA),
		convertedby2050_percent_lower = 100*sum(Converted_by2050_ha_lower)/sum(AREA_GRASS_HA),
		convertedby2050_percent_upper = 100*sum(Converted_by2050_ha_upper)/sum(AREA_GRASS_HA),
		convertedby2050_percent_private = sum(Converted_by2050_ha)/428933,
		convertedby2050_percent_lower_private = sum(Converted_by2050_ha_lower)/428933,
		convertedby2050_percent_upper_private = sum(Converted_by2050_ha_upper)/428933,
		convertedby2050_countyaverage = 100*mean(Converted_by2050_ha/AREA_GRASS_HA)
		))
	currTableB <- with(currRangeDifPerc, data.frame(
		percent_increase_on_conversion_rate = mean(percent_increase_on_background),
		percent_increase_on_conversion_lower = mean(percent_increase_on_background_lower),
		percent_increase_on_conversion_upper = mean(percent_increase_on_background_upper)
	))
	currTableC <- with(currRangeDif, data.frame(	
		additional_area_converted_by2050 = sum(Converted_by2050_ha),
		additional_area_converted_by2050_lowerCI = sum(Converted_by2050_ha_lower),
		additional_area_converted_by2050_upperCI = sum(Converted_by2050_ha_upper)
	))
	currTableD <- with(currTotalSub, data.frame(
		Area_mesic = Area_mesic_ha,
		Area_mesic_2050 = Mesic_habitat_remaining_2050_ha,
		Area_mesic_2050_lower = Mesic_habitat_remaining_2050_ha_lower,
		Area_mesic_2050_upper = Mesic_habitat_remaining_2050_ha_upper,
		Area_mesic_2050_lossperc = Mesic_habitat_loss_2050_percent,
		Area_mesic_2050_lossperc_lower = Mesic_habitat_loss_2050_percent_lower,
		Area_mesic_2050_lossperc_upper = Mesic_habitat_loss_2050_percent_upper
		))
currTableE <- with(currTotalSubnoAlf, data.frame(
		Area_mesic_noAlf = Area_mesic_ha,
		Area_mesic_2050_noAlf = Mesic_habitat_remaining_2050_ha,
		Area_mesic_2050_noAlf_lower = Mesic_habitat_remaining_2050_ha_lower,
		Area_mesic_2050_noAlf_upper = Mesic_habitat_remaining_2050_ha_upper,
		Area_mesic_2050_noAlf_lossperc = Mesic_habitat_loss_2050_percent,
		Area_mesic_2050_noAlf_lossperc_lower = Mesic_habitat_loss_2050_percent_lower,
		Area_mesic_2050_noAlf_lossperc_upper = Mesic_habitat_loss_2050_percent_upper
		))
		currTable <- data.frame(currTableA, currTableB, currTableC, currTableD, currTableE)
	return(currTable)
})
Table1 <- do.call(rbind, Table1)
Table1 <- t(Table1)
Table1 <- round(Table1, 5)
Table1<- data.frame(Table1)
names(Table1) <- c("No loss", "10% loss", "25% loss", "50% loss", "100% loss")
write.csv(Table1, "3 Model output/sage grouse impacts/Table_1_summaryofresults.csv", row.names=TRUE)

#Summary of extrapolated habitat loss 2008-2012 & 2013 to 2015
curryrs <- list("2008to2012", "2013to2015")
yrid <- list("0812", "1315")

lapply(1:2, function(i){
	actualConversionDF <- read.csv(paste0(actualConversion, sprintf("%s.csv", yrid[[i]])), header=TRUE)
	actualConversionDF <- subset(actualConversionDF, State!="CA") #drop california
	actualConversionDF <- merge(actualConversionDF, areagrassland, by="ADMIN_FIPS", all.x=TRUE)

	actualConversionTable <- with(actualConversionDF, data.frame(conversion_rate=mean(Prc),
				conversion_min=min(Prc),
				conversion_max=max(Prc), 				
				Area_Ranget0=sum(Area_Ranget0), 		
				actual_annual_conversion_ha=sum(Actual_annual_converted_ha),
				actual_annual_conversion_percent=100*sum(Actual_annual_converted_ha)/sum(AREA_GRASS_HA),  
				actual_converted_by2050_ha=sum(Actual_converted_by2050_ha),
				actual_converted_by2050_ha_percent=100*sum(Actual_converted_by2050_ha)/sum(AREA_GRASS_HA)	
				#median_mesic_ha=sum(Area_mesic_ha),
				#actual_mesic_habitat_remaining_2050_ha=sum(Actual_mesic_habitat_remaining_2050_ha, na.rm=TRUE),
				#actual_mesic_habitat_loss_2050_percent=100 - 100*sum(Actual_mesic_habitat_remaining_2050_ha, na.rm=TRUE)/sum(Area_mesic_ha)))
	))
	actualConversionTable <- t(actualConversionTable)

	write.csv(actualConversionTable, sprintf("3 Model output/sage grouse impacts/Table_1b_summaryactualconversion_%s_average.csv", curryrs[[i]]), row.names=TRUE)
})

##################
#PLOTS of MARGINAL EFFECTS
##################
#ggsave sometimes does weird things - if so open a new R and run this before doing other plots
marginplotfun <- function(dataset, currvariable) {
	p <- ggplot(dataset, aes(x=percentBLMloss, y=marginal_effect)) +
		geom_line(lwd=1)+
		geom_ribbon(aes(ymin=lowerCI, ymax=upperCI), alpha=0.1) +
		#coord_cartesian(ylim=c(0,50))+
		labs(y=sprintf("Marginal effect of %s", currvariable), x="Loss of AUM (%)") +
		theme_hc()+
		#guides(color=guide_legend(title=NULL))+
		theme(axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
		axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
		axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20))) #move title away from axis
	ggsave(filename=sprintf("4 Figures/Marginal_effect_of_%s_on_conversionprob.png", currvariable), plot=p + facet_wrap( ~ ADMIN_FIPS, ncol=8), width=8.27, height=11.69 , units="in", scale=2)	
	ggsave(filename=sprintf("4 Figures/Marginal_effect_of_%s_on_conversionprob.pdf", currvariable), plot=p + facet_wrap( ~ ADMIN_FIPS, ncol=8), device="pdf",  width=8.27, height=11.69 , units="in",scale=2)	
}

marginplotfun(rangeRentMargins, "rangeland rent")
marginplotfun(percCropMargins, "percent cropland in county")
marginplotfun(popnMargins, "county population")

#################
#FIGURE 2: Area of sagegrouse (mesic) habitat remaining by 2050 for 0-100% BLM loss
#################
#Loss curve with and without grazing restrictions. 0-100%AUM lost vs predicted area of habitat remaining (after x years), for whole study region

#as % mesic lost
p <- ggplot(totalMesicLossDF, aes(x=percentBLMloss, y=Mesic_habitat_loss_2050_percent)) +
	geom_line(lwd=1)+
	geom_ribbon(aes(ymin=Mesic_habitat_loss_2050_percent_lower, ymax=Mesic_habitat_loss_2050_percent_upper), alpha=0.1) +
	#coord_cartesian(ylim=c(0,50))+
	labs(y="Mesic habitat lost by 2050 (%)", x="Loss of AUM (%)") +
	scale_x_continuous(expand = c(0, 0.01)) + scale_y_continuous(limits=c(0,10), expand = c(0, 0.01))+ #make start at 0,0
	theme_hc(18)+
	#guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
	axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
	axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)), plot.margin = unit(c(1.5, 2, 0.5, 0.5), "lines")) #move title away from axis
ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_summedacrosslandscape.png", plot=p, width=7)	
ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_summedacrosslandscape.pdf", plot=p, width=7, device="pdf")	


#as area mesic lost
pm <- ggplot(totalMesicLossDF, aes(x=percentBLMloss, y=Mesic_habitat_loss_2050_ha/100000)) +
	geom_line(lwd=1)+
	geom_ribbon(aes(ymin=Mesic_habitat_loss_2050_ha_lower/100000, ymax=Mesic_habitat_loss_2050_ha_upper/100000), alpha=0.1) +
	#coord_cartesian(ylim=c(0,50))+
	labs(y=expression("Mesic habitat lost by 2050 (ha"~x10^{5}~")"), x="Loss of AUM (%)") +
	scale_x_continuous(expand = c(0, 0.01)) + scale_y_continuous(limits=c(0,2.5), expand = c(0, 0.01))+ #make start at 0,0
	theme_hc(18)+
	#guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
	axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
	axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)), plot.margin = unit(c(1.5, 2, 0.5, 0.5), "lines")) #move title away from axis
ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_summedacrosslandscape_ha.png", plot=pm, width=7)	
ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_summedacrosslandscape_ha.pdf", plot=pm, width=7, device="pdf")	


#as area sagebrush lost
ps <- ggplot(totalMesicLossDF, aes(x=percentBLMloss, y=Converted_by2050_ha/1000000)) +
	geom_line(lwd=1)+
	geom_ribbon(aes(ymin=Converted_by2050_ha_lower/1000000, ymax=Converted_by2050_ha_upper/1000000), alpha=0.1) +
	#coord_cartesian(ylim=c(0,50))+
	labs(y=expression("Sagebrush lost by 2050 (ha"~x10^{6}~")"), x="Loss of AUM (%)") +
	scale_x_continuous(expand = c(0, 0.01)) + scale_y_continuous(limits=c(0,2.5), expand = c(0, 0.01))+ #make start at 0,0
	theme_hc(18)+
	#guides(color=guide_legend(title=NULL))+
	theme(axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
	axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
	axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)), plot.margin = unit(c(1.5, 2, 0.5, 0.5), "lines")) #move title away from axis

ggsave(filename="4 Figures/Sagebrush_habitat_loss_byBLM_summedacrosslandscape_ha.png", plot=ps, width=7, units="in")	
ggsave(filename="4 Figures/Sagebrush_habitat_loss_byBLM_summedacrosslandscape_ha.pdf", plot=ps, width=7, units="in", device="pdf")	

#Arrange the plots together
plot1 <- arrangeGrob(ps, top=textGrob("(a)", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18)))
#fontfamily="Times Roman"
plot2 <- arrangeGrob(pm, top=textGrob("(b)", x = unit(0, "npc")
         , y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18)))
dev.off()
#grid.arrange(plot1, plot2, ncol=2) ggsave doesn't play nicely so see belwo
g <- arrangeGrob(plot1, plot2, ncol=2)
ggsave(filename="4 Figures/Sagebrush_mesic_habitat_loss_byBLM_2panel.png", plot=g, width=11, height=5, units="in")	
ggsave(filename="4 Figures/Sagebrush_mesic_habitat_loss_byBLM_2panel.pdf", plot=g, width=11, height=5, units="in", device="pdf")
ggsave(filename="4 Figures/Sagebrush_mesic_habitat_loss_byBLM_2panel.svg", plot=g, width=11, height=5, units="in")

			
#################
#FIGURE 3: Background rates of conversion
#################
#Map showing counties predicted to be most affected by conversion to cropland (Colour map just to private lands (grey out public lands))
#%change and absolute change (h)
#set up data

#make maps of ABSOLUTE area of sagebrush lost
mapfun(dataset=rangeRentvAreaAbs, currpercent=0, currfile="4 Figures/Map_absoluteconversion_%spercentBLMloss_2050hectares", fileext=".png", inset=FALSE, cornerlabel="(a) Area converted",  "Converted_by2050_ha",
			#style=c("kmeans"),
			breaks=c(-Inf, 500, 1000, 5000, 10000, Inf),
			palette=monoCols4, 
			title=c("Area converted \nby 2050 (ha)"), textNA="Data unavailable", colorNA="grey50")

			
#make maps of ABSOLUTE conversion probability

mapfun(dataset=rangeRentvAreaAbs, currpercent=0, currfile="4 Figures/Map_absoluteconversion_%spercentBLMloss_conversionprob", fileext=".png", inset=FALSE, cornerlabel="(b) Conversion rate",  "delta_conversion_prop",
			breaks=c(-Inf, 0.0005, 0.0010, 0.0015, 0.0020, Inf),
			palette=monoCols4, 
			title=c("Annual rate"), textNA="Data unavailable", colorNA="grey50")			

###SAVE HIGH RES MAPS FOR MANUSCRIPT Fig 3
#make maps of ABSOLUTE area of sagebrush lost
mapfun(dataset=rangeRentvAreaAbs, currpercent=0, currfile="4 Figures/Map_absoluteconversion_%spercentBLMloss_2050hectares", fileext=".pdf", inset=FALSE, cornerlabel="(a) Area converted",  "Converted_by2050_ha",
			#style=c("kmeans"),
			breaks=c(-Inf, 500, 1000, 5000, 10000, Inf),
			palette=monoCols4, 
			title=c("Area converted \nby 2050 (ha)"), textNA="Data unavailable", colorNA="grey50")

#make maps of ABSOLUTE conversion probability
mapfun(dataset=rangeRentvAreaAbs, currpercent=0, currfile="4 Figures/Map_absoluteconversion_%spercentBLMloss_conversionprob", fileext=".pdf", inset=FALSE, cornerlabel="(b) Conversion rate",  "delta_conversion_prop",
			breaks=c(-Inf, 0.0005, 0.0010, 0.0015, 0.0020, Inf),
			palette=monoCols4, 
			title=c("Annual rate"), textNA="Data unavailable", colorNA="grey50")				

#################
#FIGURE 4: change in conversion rate over background (counties most affected by BLM policy)
#################
#Map showing counties predicted to be most affected by conversion to cropland (Colour map just to private lands (grey out public lands))
#%change and absolute change (ha)

#make maps of additional area of sagebrush lost
mapfun(dataset=rangeRentvAreaDif, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_2050hectares", fileext=".png", inset=FALSE, cornerlabel="(b) Natural vegetation",  "Converted_by2050_ha",
			breaks=c(-Inf, 500, 1000, 5000, 10000, Inf),
			palette=seqCols2,
			title=c("Additional area\nconverted (ha)"), textNA="Data unavailable", colorNA="grey50")

#make map of percent increase on rate of sagebrush habitat loss
mapfun(dataset=rangeRentvAreaDifPerc, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_percentincreaseonrate", fileext=".png", inset=FALSE, cornerlabel="(a) Increase in conversion rate (%)", "percent_increase_on_background", 
			breaks=c(-Inf, 5, 10, 15, 25, Inf),
			palette=seqCols2, 
			title=c("Increase on \nbackground \nconversion (%)"), textNA="Data unavailable", colorNA="grey50")  

#make maps of additional area of mesic lost
mapfun(dataset=countyMesicLossDFdif, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_2050hectares_mesic", fileext=".png", inset=FALSE, cornerlabel="(c) Mesic habitat",  "Mesic_habitat_loss_2050_ha",
			breaks=c(-Inf, 50, 100, 500, 1000, Inf),
			palette=seqCols2,
			title=c("Additional area\nconverted (ha)"), textNA="No mesic", colorNA="grey50")

###SAVE HIGH RES MAPS FOR MANUSCRIPT Fig 4
#make maps of additional area of sagebrush lost
mapfun(dataset=rangeRentvAreaDif, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_2050hectares", fileext=".pdf", inset=FALSE, cornerlabel="(b) Natural vegetation",  "Converted_by2050_ha",
			breaks=c(-Inf, 500, 1000, 5000, 10000, Inf),
			palette=seqCols2,
			title=c("Additional area\nconverted (ha)"), textNA="Data unavailable", colorNA="grey50")

#make map of percent increase on rate of sagebrush habitat loss
mapfun(dataset=rangeRentvAreaDifPerc, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_percentincreaseonrate", fileext=".pdf", inset=FALSE, cornerlabel="(a) Increase in conversion rate (%)", "percent_increase_on_background", 
			breaks=c(-Inf, 5, 10, 15, 25, Inf),
			palette=seqCols2, 
			title=c("Increase on \nbackground \nconversion (%)"), textNA="Data unavailable", colorNA="grey50")  

#make maps of additional area of mesic lost
mapfun(dataset=countyMesicLossDFdif, currpercent=100, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_2050hectares_mesic", fileext=".pdf", inset=FALSE, cornerlabel="(c) Mesic habitat",  "Mesic_habitat_loss_2050_ha",
			breaks=c(-Inf, 50, 100, 500, 1000, Inf),
			palette=seqCols2,
			title=c("Additional area\nconverted (ha)"), textNA="No mesic", colorNA="grey50")			

########Other plots
			
			#make map of percent sagebrush habitat lost	
rangeRentvAreaDif$Sagebrush_habitat_lost_2050_percent <- with(rangeRentvAreaDif, delta_conversion_prop*Converted_by2050_ha)
mapfun(dataset=rangeRentvAreaDif, currpercent=50, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_percentofsagebrush", fileext=".png", inset=FALSE, cornerlabel="(a)", c("Sagebrush_habitat_lost_2050_percent"), 
			breaks=c(-Inf, 1, 2, 3, 10, 25, 100),
			labels=c("Less than 1%", "1 to 2%", "2 to 3%", "3 to 10%", "10 to 25%", "More than 25%"),
			palette=divCols, 
			title=c("Additional \narea lost \nby 2050 (%)"), textNA="Data unavailable", colorNA="grey50")  
					
#make map of percent mesic habitat lost
mapfun(dataset=countyMesicLossDFdif, currpercent=50, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_percentofmesic", fileext=".png", inset=FALSE, cornerlabel="(b)", "Mesic_habitat_lost_2050_percent", 
			breaks=c(-Inf, 1, 2, 3, 10, 25, 100),
			labels=c("Less than 1%", "1 to 2%", "2 to 3%", "3 to 10%", "10 to 25%", "More than 25%"),
			palette=divCols, 
			title=c("Additional \narea lost \nby 2050 (%)"), textNA="No mesic", colorNA="grey50")  

			#make maps of increase on conversion probability
rangeRentvAreaDif$delta_conversion_prop_perthouha <- rangeRentvAreaDif$delta_conversion_prop*1000
			
mapfun(dataset=rangeRentvAreaDif, currpercent=50, currfile="4 Figures/Map_increaseonbg_%spercentBLMloss_conversionprob", fileext=".png", inset=FALSE, cornerlabel=NA,  "delta_conversion_prop_perthouha",
			style=c("pretty"),	
			palette=divCols,
			title=c("Increase on annual \nconversion probability\n(ha/1,000 ha)"), textNA="Data unavailable", colorNA="grey50")		

#################
#FIGURE SX: MAP OF of area & percent sage grouse mesic habitat lost in each county
#################	

#make maps of area of mesic habitat lost
mapfun(dataset=countyMesicLossDF, currpercent=50, currfile="4 Figures/Map_mesic_habitat_lost_%spercentBLMloss_hectares", fileext=".png", inset=FALSE, cornerlabel="(a)", c("Mesic_habitat_loss_2050_ha"), 
			breaks=c(-Inf, 100, 500, 1000, 5000, 10000, Inf),
			palette=seqCols, 
			title=c("Area lost \nby 2050 (ha)"), textNA="No mesic", colorNA="grey50") 

#make map of percent mesic habitat lost
percentList <- list(0,50,100)
labelList <- list("(a)","(b)", "(c)")
mclapply(seq_along(percentList), function(i) {
	x=percentList[[i]]
	lab=labelList[[i]]
	mapfun(dataset=countyMesicLossDF, currpercent=x, currfile="4 Figures/Map_mesic_habitat_lost_%spercentBLMloss_percent", fileext=".png", inset=FALSE, cornerlabel=lab, c("Mesic_habitat_lost_2050_percent"), 
				breaks=c(-Inf, 1, 2, 3, 10, 25, 100),
				labels=c("Less than 1%", "1 to 2%", "2 to 3%", "3 to 10%", "10 to 25%", "More than 25%"),
				palette=divCols, 
				title=c("Area lost \nby 2050 (%)"), textNA="No mesic", colorNA="grey50")  
}, mc.cores=no_cores)

#make maps of area of mesic habitat lost #excluding conversion in alfalfa counties
mapfun(dataset=countyMesicLossDFnoAlf, currpercent=50, currfile="4 Figures/Map_mesic_habitat_lost_%spercentBLMloss_hectares_noAlfalfa", fileext=".png", inset=FALSE, cornerlabel="(b)", c("Mesic_habitat_loss_2050_ha"), 
			breaks=c(-Inf, 100, 500, 1000, 5000, 10000, Inf),
			palette=seqCols, 
			title=c("Area lost \nby 2050 (ha)"), textNA="No mesic", colorNA="grey50") 

#make map of percent mesic habitat lost #excluding conversion in alfalfa counties
percentList <- list(0,50,100)
labelList <- list("(d)","(e)", "(f)")
mclapply(seq_along(percentList), function(i) {
	x=percentList[[i]]
	lab=labelList[[i]]
mapfun(dataset=countyMesicLossDFnoAlf, currpercent=x, currfile="4 Figures/Map_mesic_habitat_lost_%spercentBLMloss_percent_noAlfalfa", fileext=".png", inset=FALSE, cornerlabel=lab, c("Mesic_habitat_lost_2050_percent"), 
			breaks=c(-Inf, 1, 2, 3, 10, 25, 100),
			labels=c("Less than 1%", "1 to 2%", "2 to 3%", "3 to 10%", "10 to 25%", "More than 25%"),
			palette=divCols, 
			title=c("Area lost \nby 2050 (%)"), textNA="No mesic", colorNA="grey50")  
}, mc.cores=no_cores)

	
#################
#FIGURE SX: SENSITIVITY ANALYSIS: COMPARE PREDICTIONS AGAINST ACTUAL (HISTORICAL) CONVERSION RATE & AREA
#################
#Plot extrapolated conversion against predicted conversion 
curryrs <- list("2008to2012", "2013to2015")
yrid <- list("0812", "1315")
numlabels <- list("(a)", "(b)")

mclapply(seq_along(yrid), function(i){
	actualConversionDF <- read.csv(paste0(actualConversion, sprintf("%s.csv", yrid[[i]])), header=TRUE)
	comparisonDF <- merge(actualConversionDF, rangeRentvAreaAbs[rangeRentvAreaAbs$percentBLMloss==0,], by="ADMIN_FIPS", all.x=TRUE)
	comparisonDF <- comparisonDF[complete.cases(comparisonDF),]
	comparisonDF$Difference_actualvpredicted_2050_ha <- comparisonDF$Converted_by2050_ha - comparisonDF$Actual_converted_by2050_ha

	print(with(comparisonDF, cor(Actual_converted_by2050_ha, Converted_by2050_ha)))
	#0.4937 2008to2012 0.5180 2013to2015
	
	p <- ggplot(comparisonDF, aes(x=Actual_converted_by2050_ha, y=Converted_by2050_ha))+
		geom_abline(slope=1, intercept=0, col="grey", lty="dashed")+	
		geom_point(aes(color=Region), size=2.5)+
		scale_colour_manual(values=divCols[c(2,3,5)]) +
		theme_classic(18)+
		#coord_cartesian(ylim=c(0,80000))+
		#guides(color=guide_legend(title=NULL))+
		theme(axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
		axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
		axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)), legend.position="bottom") +
		if(i==1){
		labs(x="Extrapolate 2008 to 2012 - sagebrush lost by 2050 (ha)", y="Model predictions - sagebrush lost by 2050 (ha)") 
		} else {labs(x="Extrapolate 2013 to 2015 - sagebrush lost by 2050 (ha)", y="Model predictions - sagebrush lost by 2050 (ha)") 
		} 
		
	labelledplot <- arrangeGrob(p, top=textGrob(numlabels[[i]], x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18)))	
	
	#ggsave(filename=sprintf("4 Figures/Scatterplot_extrapolated_vs_model_2050ha_%s.pdf", curryrs[[i]]), plot=labelledplot, device="pdf",  width=8.27, height=9.27 , units="in",scale=1)
	ggsave(filename=sprintf("4 Figures/Scatterplot_extrapolated_vs_model_2050ha_%s.png", curryrs[[i]]), plot=labelledplot,  width=8.27, height=9.27 , units="in",scale=1)

	#Make maps
	#make maps of area predicted under extrapolated values
	mapfun(dataset=comparisonDF, currpercent=curryrs[[i]], currfile="4 Figures/Map_extrapolated_hectares_%s", fileext=".png", inset=FALSE, cornerlabel=numlabels[[i]], "Actual_converted_by2050_ha",
				#style=c("kmeans"),
				breaks=c(-Inf, 1000, 2500, 5000, 10000, 20000, Inf),
				palette=seqCols, 
				title=c("Area converted (ha)"), textNA="Data unavailable", colorNA="grey50") 

	#make maps of conversion rate predicted under extrapolated values
	mapfun(dataset=comparisonDF, currpercent=curryrs[[i]], currfile="4 Figures/Map_extrapolated_conversionprob_%s", fileext=".png", inset=FALSE, cornerlabel=numlabels[[i]], "Prc",
				#style=c("pretty"),
				breaks=c(-Inf, 0.0005, 0.0010, 0.0015, 0.002, 0.003, Inf),
				labels=c("0 to 5", "5 to 10", "10 to 15", "15 to 20", "20 to 30", "More than 30"), 
				palette=seqCols, 
				title=c("Annual conversion\n(ha/10000ha)"), textNA="Data unavailable", colorNA="grey50")  

	#make difference maps of area predicted under model and extrapolated values
	mapfun(dataset=comparisonDF, currpercent=curryrs[[i]], currfile="4 Figures/Map_difference_extrapolated_vs_model_%s", fileext=".png", inset=FALSE, cornerlabel=numlabels[[i]], "Difference_actualvpredicted_2050_ha",
				breaks=c(-Inf, -5000, -500, 500, 5000, Inf),
				palette=divergingCols, 
				title=c("Difference (ha)"), textNA="Data unavailable", colorNA="grey50")
}, mc.cores=2)	
	

#################
#FIGURE SX: Map of proportion of cropland conversion on mesic habitat
#################
#make map
mapfun(dataset=countyMesicLossDF, currpercent=0, currfile="4 Figures/Map_prop_conversion_on_mesic_habitat", fileext=".png", inset=FALSE, cornerlabel=NA, "prop_conversion_on_mesic",
			style=c("pretty"),
			palette=divCols,
			title=c("Proportion of conversion \non mesic"), textNA="Data unavailable", colorNA="grey50") 
			

#################
#FIGURE SX: Sensitivity analysis of mesic lost only on non-alfalfa growing counties
#################
#Loss curve with and without grazing restrictions. 0-100%AUM lost vs predicted area of habitat remaining (after x years), for only counties where alfalfa was not the first breakout crop
p <- ggplot(totalMesicLossDF, aes(x=percentBLMloss, y=Mesic_habitat_loss_2050_percent)) +
	geom_ribbon(aes(ymin=Mesic_habitat_loss_2050_percent_lower, ymax=Mesic_habitat_loss_2050_percent_upper), alpha=0.1, fill=monoCols2[5]) +	
	geom_ribbon(data=totalMesicLossDFnoAlf, aes(ymin=Mesic_habitat_loss_2050_percent_lower, ymax=Mesic_habitat_loss_2050_percent_upper), alpha=0.1, fill="darkblue") +	
	geom_line(lwd=1, col=monoCols2[5])+
	geom_line(data=totalMesicLossDFnoAlf, aes(x=percentBLMloss, y=Mesic_habitat_loss_2050_percent), lwd=1, col="darkblue") +
	#coord_cartesian(ylim=c(0,50))+
	labs(y="Mesic habitat lost by 2050 (%)", x="Loss of AUM (%)") +
	scale_x_continuous(expand = c(0, 0.01)) + scale_y_continuous(limits=c(0,7.5), expand = c(0, 0.01))+ #make start at 0,0
	theme_hc(18)+
	theme(legend.position="bottom", axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), 
	axis.ticks.x = element_line(colour="grey"), axis.ticks.y = element_line(colour="grey"), #change tick colour
	axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)), plot.margin = unit(c(1.5, 2, 0.5, 0.5), "lines")) #move title away from axis
ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_noalfalfacounties.png", plot=p, width=7)	
#ggsave(filename="4 Figures/Mesic_habitat_loss_byBLM_noalfalfacounties.pdf", plot=p, width=7, device="pdf")		

#make maps of area of mesic habitat remaining #excluding conversion in alfalfa counties
mapfun(dataset=countyMesicLossDFnoAlf, currpercent=50, currfile="4 Figures/Map_mesic_habitat_remaining_%spercentBLMloss_hectares_noAlfalfa", fileext=".png", inset=FALSE, cornerlabel="(a)", "Mesic_habitat_remaining_2050_ha", 
			breaks=c(-Inf, 100, 500, 1000, 5000, 10000, Inf),
			palette=rev(divCols), 
			title=c("Mesic remaining \nby 2050 (ha)"), textNA="No mesic", colorNA="grey50") 

#make map of percent mesic habitat remaining #excluding conversion in alfalfa counties
mapfun(dataset=countyMesicLossDFnoAlf, currpercent=50, currfile="4 Figures/Map_mesic_habitat_remaining_%spercentBLMloss_percent_noAlfalfa", fileext=".png", inset=FALSE, cornerlabel="(b)", "Mesic_habitat_remaining_2050_percent", 
			style=c("kmeans"),
			palette=rev(divCols), 
			title=c("Mesic remaining \nby 2050 (%)"), textNA="No mesic", colorNA="grey50")  



###################
#Map of first breakout crop in each county
firstcrop <- read.csv("1 Inputs/Cropland and pasture/First_crop_byCounty_2008to2012_aspercent.csv", header=TRUE)
firstcrop$most_common_breakout_crop[firstcrop$most_common_breakout_crop=="Spring.Wheat"] <- "Wheat"
firstcrop <- droplevels(firstcrop)
cropShp <- merge(sgCountiesShp, firstcrop, by.x='fips', by.y='ADMIN_FIPS') 


	cropmap <- tm_shape(cropShp) +
		tm_fill("most_common_breakout_crop", palette="Pastel2", textNA="No crop conversion", colorNA="grey95", title="Most common \nbreakout crop") +
		tm_style_white()+
		tm_layout(title=NA, title.size=1.6, inner.margins=c(.1,.12, .02, .1), outer.margins=c(0,0,0,0), frame=TRUE,
				legend.title.size = 1,legend.text.size = 0.6, legend.position = c("right","bottom"), legend.bg.color = "white", legend.bg.alpha = 1, asp=1.25) +
		#tm_shape(publicShp) + 
		#	tm_fill(col="white")+
		tm_shape(cropShp) +
			tm_borders(lwd=0.5, col = "grey60", alpha=0.5) + #county boundaries			
		tm_shape(usStatesShp) +
			tm_borders(lwd=1, col = "grey70", alpha = 0.5)
			#tm_text("STATE", col="grey50")

save_tmap(cropmap , filename="4 Figures/Map_mainbreakoutcrop_incounty.png", outer.margins=0.01, width=7)

#################
#FIGURE 4: TIME TO LOSS of LEKS of (hockeystick) sage grouse habitat in each county
#################
#Joe Smith found that 96% of active leks are located in landscapes with less than 15% cropland
#For each ten percentage point increase in cropland, an associated 54% decrease in lek density.
#time to loss of sg leks from county is when %cropland in county ==15%
#area of sage grouse habitat affected is area of sage grouse habitat in LCC1to6 in that county
#% sage grouse habitat lost = area sg habitat in LCC1to6/(total area sg habitat in county)

#OR
#sg habitat lost = lek density in LCC1to6 * hockey stick fn(% converted by 2050) / (total area sg habitat in county)


