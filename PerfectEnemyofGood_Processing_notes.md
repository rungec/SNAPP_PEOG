# Processing notes for Perfect enemy of good paper
SNAPP Better Land Use - Sage Grouse project

## Rental rates datasets investigated
### USDA NASS land rents:
raw data was processed and summarised using script *PerfectEnemyofGood_model_prelimtests.r*  

*Source:* USDA Survey. Data from https://quickstats.nass.usda.gov/results/8275FFAE-5319-3417-AA97-F95A3384A3AD 
The Cash Rents Survey provides the basis for county estimates of the cash rent paid for irrigated cropland, non-irrigated cropland, and pasture.  The 2008 Farm Bill mandated that NASS provide mean rental rates for all counties with 20,000 acres of cropland plus pasture.
The county level Cash Rents Survey is conducted bi-annually in all states, except Alaska, beginning with a release in September, 2016. U.S. and state estimates are released in August every year using the June Area Survey. Previously, an annual survey was conducted from 2008 until 2014.  All qualifying counties in these states are represented in the sample.  The target population is all farms and ranches that have historically rented land on a cash basis for any of the three land use categories.  Land rented for a share of the crop, on a fee per head, per pound of gain, by animal unit month (AUM), rented free of charge, or land that includes buildings such as barns are excluded from the survey. More info: 
https://www.nass.usda.gov/Surveys/Guide_to_NASS_Surveys/Cash_Rents_by_County/
Approximately 240,000 farms and ranches across the United States are contacted for their total acres operated and acres rented for cash for each land use category (irrigated cropland, non-irrigated cropland, and permanent pasture) for the current year.  For each land use category with positive acres, respondents are given the option of reporting rent per acre or total dollars paid.
Uses
The Farm Service Agency (FSA) uses cash rent county estimates to determine market-based rates in administering USDA programs, such as the Conservation Reserve Program (CRP).  Other state and federal government agencies, universities, and research organizations use these data for other forms of economic analysis.  The data provide farmers and ranchers with current information about rental rates in their county and are available for their use in making decisions regarding renting and leasing farmland.
Frequency
The Cash Rents Survey is conducted annually beginning with the first mail out in mid-February.
Methods
Data collection for the Cash Rents Survey is conducted by mail and telephone.  The initial mail-out occurs in mid-February followed by a second mailing in mid-March.  Non-response phone follow-up is conducted from April to July from NASS’ Data Collection Centers.  Some field follow up is conducted to insure adequate response rates for the survey.  Over 35 percent of the survey is collected by mail.
Data are summarized to provide the mean cash rental rate for each land use category in each district and county.  District and county mean rates must reconcile to a previously published state mean cash rental rate for each category.  State level mean cash rents are estimated from the Cash Rents Survey.  Samples for this survey are drawn with a county-level stratified design.  An additional indication is also available from the June Area Survey.  This is a probability sample of land segments selected from the complete NASS area sampling frame stratified and sampled by intensity of agriculture.

*Folder:* 1 Inputs/Rental Rates/USDA NASS land rent by county/  

*Categories:* 
RENT, CASH, PASTURELAND - EXPENSE, MEASURED IN $ / ACRE  
RENT, CASH, CROPLAND, NON-IRRIGATED - EXPENSE, MEASURED IN $ / ACRE  
RENT, CASH, CROPLAND, NON-IRRIGATED - EXPENSE, MEASURED IN $ / ACRE  
>+ *Farm:* Any establishment from which $1,000 or more of agricultural products were sold or would normally be sold
during the year. Government payments are included in sales.
Farm real estate value: The value at which all land and buildings used for agriculture production including dwellings,
could be sold under current market conditions, if allowed to remain on the market for a reasonable amount of time.
>+ *Cropland value:* The value of land used to grow field crops, vegetables or land harvested for hay. Land that switches
back and forth between cropland and pasture should be valued as cropland. Hay land, idle cropland and cropland enrolled
in government conservation programs should be valued as cropland.
>+ *Irrigated cropland value:* The value of land that normally receives or has the potential to receive water by artificial
means to supplement natural rainfall. Irrigated cropland may consist of both land that will or will not be irrigated during
the current year, but still has the facilities and equipment to do so. Irrigation facilities and equipment such as wells,
pumps, canals, ditches, reservoirs, lakes, tanks, ponds, rivers, streams or creeks are usually present or on nearby acres.
>+ *Non-irrigated cropland value:* The value of land that only receives water by natural rainfall.
>+ *Pasture, grazing and grassland value:* The value of land that is normally grazed by livestock. Pasture does not need to
have livestock grazing on it at the time of interview or during the current year in order to be valued as pasture or grazing
land.

*Years available:* rent data is available annually for 2008 onwards, though 2008 data is sparse. ~730 counties have net returns data for rangeland & cropland for 2008, ~2500 counties have net returns data for 2009-2012.  

*Issues with dataset:* Classification into cropland & pasture are not comparible with remotely sensed land use datasets. It is not clear whether 'pasture' (as in land used to grow hay for grazing) is classified in this dataset as 'Pasture' or if it is included in the 'Cropland' classification. The 'Pasture' class of this dataset signifies rangeland.  

*Summary statistics* of dataset are available in:  
> 1_PEOG_grazing_economics\4 Figures

*CPI-adjustment*
Values were adjusted for average annual CPI CPI-All Urban Consumers (Current Series) CUUR0000SA0 from US Department of Labor Bureau of Labor Statistics (https://www.bls.gov/cpi/). Correction factors can be found here:  
> CPI_adjusted.csv  

###Jianhong Mu net returns for cropland, pasture & rangeland
raw data was processed and summarised using script *PerfectEnemyofGood_model_prelimtests.r*  

*Source:* An Empirical Analysis of Climate Uncertainty and Land-use Transitions in the U.S. Pacific and Mountain Regions, Jianhong E. Mu, Christopher Mihiar, David J. Lewis, Benjamin Sleeter, and John T. Abatzoglou, Selected Paper prepared for presentation for the 2016 Agricultural & Applied Economics Association, Boston, MA, July 31-August 2.  

*Folder:* 1 Inputs/Rental Rates/Net returns from Jianhong Mu/  

*Categories:* Cropland, Pasture, Rangeland  

*Years available:* 1978-present, annually (see note below).  ~990 counties have net returns data for rangeland & cropland across the time period 2007-2011.  

*Issues with dataset:* 	Two issues  
1. In many counties in our study region the net returns for cropland are quite negative (and less than the net returns for rangeland or pasture)  
2. Economic data (farm expenses) has been portioned out into crop,range,pasture using NRI land use transitions (which are not the best source of data, summary above), and assumes land use remains static between 2007-2011 This might not be such a big issue, the economic data is portioned out according to the proportional area of each land use type in the landscape. The proportion of land in any of these land classes that is changing in a given year is small (typically <1%), so innacurracy arising from this (in the net return data) should be correspondingly small.  

*Notes:* Net returns were calculated from annual county-level data of farm income and expenses for crop and livestock production collected from the Bureau of Economic Analysis (BEA). Total farm production revenue includes sales from crop production, sales from livestock production and total governmental payments and minus payment for conservation reserve programs. The per acre cropland net return is adjusted for the acreas of cropland in a county. Because sales of livestock products were not spilt by pasture- or rangeland, they did additional calculations to get approximated net return per acre for pasture- and rangeland. Basically, the per acre sales of livestock products are portioned out according to the area of rangeland & pasture in a county, using either NRI or Lark estimates of areas of each land use in a county. Further details can be found in the reference above.  

*Summary statistics* on comparisons between the data adjusted using NRI vs the data adjusted using Lark are available in:  
> Comparison_of_net_returns_adjusted_using_NRI_vs_using_Lark.txt
> Net returns_differencebetweenLarkNRI_byYear.png
I concluded that the NRI dataset is better, as JMu has done some processing to reign in outliers in the dataset, and because I think the fact that the pasture areas from Donnelly only are in sage grouse habitat is skewing the ratios of rangeland:pastureland in the Lark-adjusted data.  

## Cropland conversion rates datasets investigated
Conversion Probabilty is the net change from (pasture or cropland) to rangeland in a given county & year, divided by the area of rangeland in that county & year.  

### Remotely sensed conversion rates
Using *PerfectEnemyofGood_prelimmodels.r*  
Combined annual area of cropland, pasture, rangeland from data sourced from Tyler Lark & Patrick Donnelly into a single dataset (areas are in hectares) and calculated conversion rates.  
Areas are calculated on only private lands that are not forest, water, urban or barren in the 2015 CDL, for all Land capability classes. Areas are in hectares. See *ACR-grassland-conversion/ProcessingNotes.md* for details on pre-processing of the Lark dataset.  
> CombinedCropRangePasture_areas_LarkDonnelly_and_conversionprobabilities_Allcounties.csv
> CombinedCropRangePasture_areas_LarkDonnelly_and_conversionprobabilities_sgcounties.csv

Input datasets are:  
#### Lark CDL
*Source:* Tyler Lark lark@wisc.edu  
[check whether this is the preferred citation]. Lark, T.J., Salmon, J.M. & Gibbs, H.K. (2015). Cropland expansion outpaces agricultural and biofuel policies in the United States. Environmental Research Letters, 10, 44003.  
*Categories:* Cropland, Rangeland - follows Cropland Data Layer (CDL) categories  
*Years available:* 2008-2015 
*Extent:* Cont. USA  
*Summary* available in reference, and in NCEAS_Postdoc\P4 ACR revised methods\Analysis\figures  
*Notes:* This dataset is an improvement on the CDL that addresses a crop area underestimation bias within the CDL that needs to be corrected when calculating changes to crop and cropland areas.  

Note: Classification is difficult between grassland and alfala in wetter areas. In this dataset, something classed as grassland initially that converts to alfalfa is classed as NONCROP. grassland that converts to alfalfa, then to another crop (corn, wheat) is classed as CROP. land that is alfalfa continuously is classed as STABLE CROP.

Layers (2008 to 2012, 56m resolution) are:
>+ *Multitemporal_Results_FF2.tif* 1=stable noncrop, 2= stable crop, 3= converted to crop, 4= abandoned, 5=intermittent cropland
>+ *year_from_crop_ff2.tif* value is year cropland was abandoned
>+ *year_to_crop_ff2.tif* value is year land was converted to cropland
>+ *class_before_crop.tif* value corresponds to the CDL classification for the last land use before land was converted to cropland
>+ *crop_before_noncrop.tif* value corresponds to the CDL classification for the first crop grown on a pixel after it was converted to cropland  
>+ *first_crop_class.tif* value corresponds to the CDL classification for the first crop grown on a pixel after it was converted to cropland  
>+ *first_noncrop_class.tif* value corresponds to the CDL classification for the first land use after cropland was abandoned

Layers (2013 to 2015, 30m resolution) are stored in xp_update.gdb:
>+ *mtc* (conversion type) 
>+ *ytc* (year to crop)
>+ *ytc_bfc* (crop before conversion)
>+ *ytc_fc* (first crop)

*Data processing:*  
Areas that are forest, developed or water were masked from the Lark CDL rasters, based on the 2015 CDL. These rasters were then summarised (tabulate by area) using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to Multitemporal_Results_FF2.tif to give tables of how much land was converted or present as cropland in each year, and this analysis is detailed in:  
https://github.nceas.ucsb.edu/rungec/ACR-grassland-conversion/ProcessingNotes.md  

Conversion probabilities are calculated for private (non-government) land that is not classified as water, forest or developed in the 2015 CDL (i.e. "grassland", includes only pixels classified as shrubland, wetlands, native grassland and pasture/hay).  

The 2012 to 2015 data from Tyler Lark was at 30m resolution, 2008 to 2012 data was at 56m resolution. The 2013-2015 data has false mapping of alfalfa/fallow (notes from Tyler: *We found 2 counties with substantial false change identified, which you might want to exclude from your analysis: 
18087  Lagrange, IN
48249  Jim Wells, TX
In these spots, the CDL shows lots of false conversion to alfalfa (in Lagrange) or to Fallow (in Jim Wells).   There could also be a few other outlier counties with similar issues but less pronounced. What I'd recommend: use the "first_crop_class" layer to exclude or mask out all areas of identified conversion where the first crop class was Alfalfa (36) or Fallow (61).   That should fix both of those problem counties, as well as any others with the same issue.  We had to do something similar to our old 2008-2012 data, so shouldn't affect compatibility at all.*) Using con in ArcGIS, I masked out areas where the first crop class was alfalfa or fallow.     

Tables used (which can be found in P4 ACR revised methods\Analysis\tables\land use by area) were:  
> LarkCDL_GrasslandPrivateArea_byCounty.csv
> LarkCDL_GrasslandPrivateYearFromCrop_byCounty.csv
> LarkCDL_GrasslandPrivateYearToCrop_byCounty.csv

##### Mask Forest, Water and Developed land from Lark's CDL data
*manual in ArcGIS*  
Conversion probabilities are calculated for land that is not classified as water, forest or developed in the 2015 CDL (i.e. grassland mask includes only pixels classified as shrubland, wetlands, native grassland and pasture/hay).  

Use 'Con' tool (conditional). If Class_Name in 2015 CDL layer 2015_30m_cdls.img is
"Class_Name" = 'Clouds/No Data' OR "Class_Name" = 'Deciduous Forest' OR "Class_Name" ='Developed'OR "Class_Name" = 'Developed/High Intensity'OR "Class_Name" = 'Developed/Low Intensity' OR "Class_Name" ='Developed/Med Intensity'OR "Class_Name" = 'Developed/Open Space' OR "Class_Name" ='Evergreen Forest' OR "Class_Name" ='Forest'OR "Class_Name" = 'Mixed Forest' OR "Class_Name" ='Open Water'OR "Class_Name" = 'Perennial Ice/Snow' OR "Class_Name" ='Water'
then reclassify Multitemporal_Results_FF2.tif/mtc as value 15 (no data). Snap to and extent of Multitemporal_Results_FF2.tif/mtc
> LarkCDL_grassland.tif
> LarkCDL_grassland_2012_2015.tif

Use 'Con' tool (conditional). If Class_Name in 2015 CDL layer 2015_30m_cdls.img is water, forest, or developed
then reclassify year_from_crop_ff2.tif/ytc as "" no data. Snap to and extent of year_from_crop_ff2.tif
> LarkCDL_yearfromcrop_grassland.tif

Use 'Con' tool (conditional). If Class_Name in 2015 CDL layer 2015_30m_cdls.img is water, forest, or developed then reclassify year_to_crop_ff2.tif/ytc as value 65535 (or no data). Snap to and extent of year_to_crop_ff2.tif/ytc  
> LarkCDL_yeartocrop_grassland.tif
> LarkCDL_yeartocrop_grassland_2012_2015.tif

##### Drop areas where first crop was alfalfa or fallow from 2013_2015 data
*manual in ArcGIS* 
Use 'Con' tool (conditional). If Value in ytc_fc is 36 or 61 (alfalfa/fallow) then reclassify LarkCDL_yeartocrop_grassland_2012_2015.tif/LarkCDL_grassland_2012_2015.tif as no data. Snap to and extent of LarkCDL_grassland_2012_2015.tif
> LarkCDL_grassland_2012_2015_noalfalfafallow.tif
> LarkCDL_yeartocrop_grassland_2012_2015_noalfalfafallow.tif

##### Make rasters of conversion for each of the 8 land capability classes
Reclassify all areas outside a given capability class as no data.  
*manual in ArcGIS*  
Use 'Con' tool (conditional). If Class_Name in Land Capability Class layer LCC_100m.tif (for 2008-2012 data)/LCC_100m_albersNAD83.tif (for 2013-2015 data) is "Value" <> 1 (or "Value" >= 6 etc) then reclassify LarkCDL_grassland.tif/LarkCDL_grassland_2012_2015_noalfalfafallow.tif as value 16 (or no data). Snap to and extent of Multitemporal_Results_FF2.tif/LarkCDL_grassland_2012_2015_noalfalfafallow.tif etc  
> LarkCDL_grassland_LCC1.tif
> LarkCDL_grassland_LCC2.tif
> LarkCDL_grassland_LCC3.tif
> LarkCDL_grassland_LCC4.tif
> LarkCDL_grassland_LCC5.tif
> LarkCDL_grassland_LCC6.tif
> LarkCDL_grassland_LCC7.tif
> LarkCDL_grassland_LCC8.tif
> LarkCDL_grassland_LCC1to4.tif
> LarkCDL_grassland_LCC1to6.tif
> LarkCDL_grassland_LCC7or8.tif
> LarkCDL_grassland_LCC5or6.tif
> LarkCDL_grassland_2012_2015_noalfalfafallow_LCC1to4.tif
> LarkCDL_grassland_2012_2015_noalfalfafallow_LCC1to6.tif
> LarkCDL_grassland_2012_2015_noalfalfafallow_LCC5or6.tif
> LarkCDL_yearfromcrop_grassland_LCC1to4.tif
> LarkCDL_yearfromcrop_grassland_LCC1to6.tif
> LarkCDL_yearfromcrop_grassland_LCC7or8.tif
> LarkCDL_yeartocrop_grassland_LCC1to4.tif
> LarkCDL_yeartocrop_grassland_LCC1to6.tif
> LarkCDL_yeartocrop_grassland_LCC7or8.tif
> LarkCDL_yeartocrop_grassland_LCC5or6.tif
> LarkCDL_yeartocrop_grassland_2012_2015_noalfalfafallow_LCC1to6.tif
> LarkCDL_yeartocrop_grassland_2012_2015_noalfalfafallow_LCC1to4.tif
> LarkCDL_yeartocrop_grassland_2012_2015_noalfalfafallow_LCC5or6.tif

##### Calculate area of grassland converted to cropland in each county
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to Multitemporal_Results_FF2.tif  
Classes of LarkCDL_grassland.tif are 1=stable noncrop, 2= stable crop, 3= converted to crop, 4= abandoned, 5=intermittent cropland , 15=forest,water or developed  
Joined to PADUSCBIv2_Private_Land_only_AllCounty_diss by ADMIN_FIPS and exported table    
> input: LarkCDL_grassland.tif / LarkCDL_grassland__2012_2015_noalfalfafallow.tif
> output:LarkCDL_GrasslandPrivateArea_byCounty.csv / LarkCDL_GrasslandPrivateArea_byCounty_2012_2015.csv 
Areas are in m2. The sum of values 0 to 5 is the sum area of private grassland or cropland in each county.

##### Calculate area of grassland converted to cropland in each county in each year
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of  using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to Multitemporal_Results_FF2.tif  
Joined to PADUSCBIv2_Private_Land_only_AllCounty_diss by ADMIN_FIPS and exported table  
> input: LarkCDL_yeartocrop_grassland.tif / LarkCDL_yeartocrop_grassland_2012_2015_noalfalfafallow.tif
> output: LarkCDL_GrasslandPrivateYearToCrop_byCounty.csv  / LarkCDL_GrasslandPrivateYearToCrop_byCounty_2012_2015.csv
Areas are in m2 of private grassland converted to cropland in each year in each county.

##### Calculate area of cropland converted to other use in each county in each year
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of LarkCDL_yearfromcrop_grassland.tif using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to Multitemporal_Results_FF2.tif  
Joined to PADUSCBIv2_Private_Land_only_AllCounty_diss by ADMIN_FIPS and exported table 
> input: LarkCDL_yearfromcrop_grassland.tif
> output: LarkCDL_GrasslandPrivateYearFromCrop_byCounty.csv  
Areas are in m2 of private grassland converted to cropland in each year in each county.

##### Calculate area of each land capability class in each county (all land)
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of LCC_100.tif using USGS_County_boundaries_USContinental_albers.shp ADMIN_FIPS as zones, and with cell size set to 100m 
Classes of LarkCDL_grassland.tif are 1 to 8 (1-6 suitable, 7-8 unsuitable)
Joined to USGS_County_boundaries_USContinental_albers.shp by ADMIN_FIPS and exported table    
> Areaof_LandCapability_byCounty.csv  

##### Calculate area of each land capability class in each county (private unprotected land only)
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of LCC_100.tif using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to 100m 
Classes of LarkCDL_grassland.tif are 1 to 8 (1-6 suitable, 7-8 unsuitable)
Joined to PADUSCBIv2_Private_Land_only_AllCounty_diss by ADMIN_FIPS and exported table    
> Areaof_LandCapability_byCounty_PrivateLandOnly.csv  

#### Calculate area of grassland in each land capability class converted to cropland in each county
*manual in ArcGIS*  
Tabulate by area (Zonal toolbox) of LarkCDL_grassland_LCC1.tif (etc) using PADUSCBIv2_Private_Land_only_AllCounty_diss.shp ADMIN_FIPS as zones, and with cell size set to LarkCDL_grassland_LCC1.tif (etc) 
Classes of LarkCDL_grassland_LCC1.tif are 1=stable noncrop, 2= stable crop, 3= converted to crop, 4= abandoned, 5=intermittent cropland , 15=forest,water or developed, 16=other LCC  
Joined to PADUSCBIv2_Private_Land_only_AllCounty_diss by ADMIN_FIPS and exported table    
> LarkCDL_GrasslandPrivateArea_LCC1_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC2_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC3_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC4_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC5_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC6_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC7_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC8_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC1to4_byCounty.csv  
> LarkCDL_GrasslandPrivateYearFromCrop_LCC1to4_byCounty.csv  
> LarkCDL_GrasslandPrivateYearToCrop_LCC1to4_byCounty.csv 
> LarkCDL_GrasslandPrivateArea_LCC1to6_byCounty.csv  
> LarkCDL_GrasslandPrivateYearFromCrop_LCC1to6_byCounty.csv  
> LarkCDL_GrasslandPrivateYearToCrop_LCC1to6_byCounty.csv  
> LarkCDL_GrasslandPrivateArea_LCC7or8_byCounty.csv 
> LarkCDL_GrasslandPrivateYearFromCrop_LCC7or8_byCounty.csv 
> LarkCDL_GrasslandPrivateYearToCrop_LCC7or8_byCounty.csv 
> LarkCDL_GrasslandPrivateYearToCrop_LCC1to4_byCounty_2012_2015.csv  
> LarkCDL_GrasslandPrivateYearToCrop_LCC1to6_byCounty_2012_2015.csv  
> LarkCDL_GrasslandPrivateYearToCrop_LCC5or6_byCounty_2012_2015.csv  
> LarkCDL_GrasslandPrivateArea_LCC1to4_byCounty_2012_2015.csv  
> LarkCDL_GrasslandPrivateArea_LCC1to6_byCounty_2012_2015.csv  
> LarkCDL_GrasslandPrivateArea_LCC5or6_byCounty_2012_2015.csv  
Areas are in m2. The sum of values 0 to 5 is the sum area of private grassland or cropland in each LCC in each county.

I then used these tables made for the ACR project to get the area of cropland in each county in each year, for use in the PEOG project. Calculations were performed in Cropland_and_Rangeland_Area_byYear_ha.xlsx & LarkCDL_GrasslandPrivateArea_byLCC_byCounty_andYear.xlsx

The area of cropland in each year is the area of cropland at the beginning of the time period (2008; the sum of the stable cropland, plus any cropland added in that year and previous years.  

The area of rangeland in each year is the area of rangeland at the beginning of the time period (2008; the sum of stable non-crop, plus the area of cropland that was lost between 2008 and 2012), minus any rangeland converted to cropland in that year and previous years. 

These values are for private, unprotected land, and excludes land that is classified in the 2015 CDL as water, forest or developed. ie. only represents cropland transitions to or from rangeland.  

We ignore cropland that is abandoned or intermittent, though in 1796 counties the area of intermittent cropland is greater than the area of land converted to cropland (vs 1298 smaller); in 1381 counties the area of land abandoned is greater than the area of land converted to cropland (vs 1722 smaller). Even with investment in restoration is unlikely at this land would be restored to sage grouse habitat within the time periods we consider.

(1) it is unlikely to become sage grouse habitat in that time scale we are looking at. (2) Abandoned land may be used for pasture, or may simply be left fallow.

*no longer done*{Abandoned cropland was assumed to be converted to rangeland instantaneously (ie without requiring reveg, and it was assumed that all land lost from cropland goes to rangeland. This assumption is reasonable, given that >90% of counties more than 50% of abandoned cropland is converted to rangeland rather than developed or other. See Class_after_crop_bycounty.xlsx for analysis.}  

input: Cropland_and_Rangeland_Area_byYear_ha.xlsx
output: (hectares)
> Rangeland_Area_byYear_ha.csv
> Cropland_Area_byYear_ha.csv

input: LarkCDL_GrasslandPrivateArea_byLCC_byCounty_andYear.xlsx
output: (hectares)
> Rangeland_Area_byYear_LCC1to4.csv
> Rangeland_Area_byYear_LCC1to6.csv
> Rangeland_Area_byYear_LCC7or8.csv
> CroplandArea_byYear_LCC1to4.csv
> CroplandArea_byYear_LCC1to6.csv
> CroplandArea_byYear_LCC7or8.csv



#### Summarise & combine LULC conversion rates
*PerfectEnemyofGood_mergeData.r*
Combined conversion probabilites for all land, LCC1to4, LCC1to6, LCC7or8 into a combined dataset for the purposes of PEOG models. 
> Allcounties_LarkConversion_byLCC_GrasslandPrivateArea.csv


#### Donnelly mesic sage grouse habitat
*Source:* Patrick Donnelly, [citation to follow, still in prep]
*Years available:* 1984-present  
*Extent:* Western US, areas overlapping sage grouse habitat.  
*Notes:* This is a remotely sensed dataset of temporal changes in mesic sage grouse habitat.  

Filters he applied to the County_mesic.csv dataset:  
county  
private ownership  
wet meadows (i.e. hay pastures)  
year  
viable late brood rearing habitat (ha)  
The areas represent hay pastures that are <200m from sagebrush

In total, there is 6410780 acres (average) and 639633 acres (median) of sage brush habitat in wet-meadows on private land in all states (not just in counties we looked at) (from County_mesic.csv). There is 3,883,000ha of any type of mesic habitat (alfalfa, mesic rangeland, riparian, wet meadow) on private land in the study region (from MESIC_private_Lark_diss.shp below), compared to 48,230,000 ha private land. 5,697,680ha is cropland or converted to cropland. 219,885 ha is converted to cropland. 23,512 ha of that was on mesic habitat.

Processing:
*Manual in ArcGIS*
First merged these four regional .shps (regions are Colombia Basin, Great Plains, Great Basin, Rocky Mountains), 
then union with private lands SGCounties_private.shp, followed by clip to SGCounties_private.shp.  
I then clipped Multitemporal_Results_FF2.tif (land use, from 2008 to 2012 Lark CDL data) to study region, then converted it to a polygon, reprojecting Albers equal area (USGS version) to Transverse Mercator (to match mesic shps) (output saved to temp file) then intersect the polygon version of Multitemporal_Results_FF2.tif and MESIC_Merge_private. Finally I dissolved by ADMIN_FIPS, TYPE (of mesic habitat) & grid_code (the value from Multitemporal_Results_FF2.tif).

> 1_PEOG_grazing_economics\1 Inputs\Sage grouse\Mesic_habitat_Donnelly\
> Mesic_Claire.gdb/MESIC_Merge.shp
> Mesic_Claire.gdb/MESIC_Merge_private.shp
> MESIC_private_Lark.shp
> MESIC_private_Lark_diss.shp
> MESIC_private_Lark_diss.csv

This .csvs was then formatted using script *PerfectEnemyofGood_createMesicInputs.r*
> MESIC_private_Lark_byCounty_wide.csv


### Survey-based conversion rates
#### NRI
Land use for each plot in the NRI survey is reported every 5 years (...2002, 2007, 2012). Data is obtained from a subset of farms across the US. Because it is not a full sample, and the probability of any given farm converting use is low, the conversion rates estimated from this dataset are subject to a high uncertainties. Though there are issues with distinguishing similar crop types by remote sensing, remotely sensed data would be expected to give a more accurate estimate of cropland conversion rates, when aggregated to a generic 'cropland' class. 

#### USDA
I havent looked at this data in detail. It would be subject to similar issues to NRI.

## Other datasets
### Land ownership
*Source:*  
> https://catalog.data.gov/dataset/usgs-small-scale-dataset-1-1000000-scale-county-boundaries-of-the-united-states-201403-filegdb  
> http://consbio.org/products/projects/pad-us-cbi-edition  

*Notes:* Private land: PADUS_CBIEdition_V2 selected by attributes only private lands, merged with county boundaries from USGS_digital_boundaries, then dissolved by county.     
> PADUSCBIv2_Private_Land_only_AllCounty_diss.shp  
Erase strict protected areas, dissolved then intersect with county boundaries. Strict protected areas defined as either GAP_Sts 1 or 2, or IUCN category I-VI  
> PADUSCBIv2_11sgStates_nostrictPAs_byCounty.shp  

### Land use
National Land Cover Database 2011 (NLCD 2011) is the most recent national land cover product created by the Multi-Resolution Land Characteristics (MRLC) Consortium. NLCD 2011 provides - for the first time - the capability to assess wall-to-wall, spatially explicit, national land cover changes and trends across the United States from 2001 to 2011. As with two previous NLCD land cover products NLCD 2011 keeps the same 16-class land cover classification scheme that has been applied consistently across the United States at a spatial resolution of 30 meters. NLCD 2011 is based primarily on a decision-tree classification of circa 2011 Landsat satellite data.  
Preferred NLCD 2011 citation: Homer, C.G., Dewitz, J.A., Yang, L., Jin, S., Danielson, P., Xian, G., Coulston, J., Herold, N.D., Wickham, J.D., and Megown, K., 2015, Completion of the 2011 National Land Cover Database for the conterminous United States-Representing a decade of land cover change information. Photogrammetric Engineering and Remote Sensing, v. 81, no. 5, p. 345-354  
Accessed 10 December 2015  
> http://www.mrlc.gov/nlcd2011.php  

### AUM data
Rangeland productivity was obtained from gSSURGO. The gSSURGO original data downloaded from the USDA website are in ESRI's FileGeodatabase, we had to export them to outside de FGDB to perform the tasks in R. We exported the Soil's rasters to Geotiff (.tif) and the atribute tables to DBF. 
Using scripts *JoinMosaic_SGMZ.R, JoinMosaic_SGMZ_gdbExtractor.py*  
This script reads the *soil raster* and the *'component'* table for each state as inputs, selects only the favourable, normal and unfavourable forage productivity (prodfav, prodnorm, produnfav) field, performs the join of both data based on the 'mukey' field, creates a raster for each state with values of productivity and (optionally) mosaic all states into a larger raster("SGMZ_Mosaic.tif").  

This data is summarised by county using *PerfectEnemyofGood_summariseForageAUM.r*
Summary by county, by ownership and whether or not it falls in sage grouse habitat (overlaid with SGCounties_sgPACsBr_rastint.tif)
> AUM stocking rates\forageZonal_SGCountiesAllLand

Summary by county, by ownership and LCC, exclusing forest, developed & urban, and strict protected areas
(overlaid forage rasters by PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp and a masked lcc raster, based on Id field. Values in Id: first digit=LCC class; second digit=private(1) or public (0); last 5 digits=ADMIN_FIPS)
> AUM stocking rates\forageZonal_11states_grassland_byLCC\

These forage levels were then converted to AUM: using *PerfectEnemyofGood_summariseForageAUM.r*
AUM is defined as  (see AUM_calculation_notes.xlsx for more details)
> AUM_peracre = round(foragedf[,x]*0.25/(30*30.5), 3)   
> *Carrying capacity in AUM, grazin efficiency 25%, 30lb air dried feed per day, 30.5days per month* 
> AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_raw.csv

Summarised for AllLCC, LCC1to6, LCC1to4 and formatted to match datasets used in econometric models: using *PerfectEnemyofGood_summariseForageAUM.r*
> AUM_and_forage_estimates_fromgSSURGO_byStateCountyLCCPublicPrivate_final.csv
 

*Key reference:* Butler LD, Cropper JB, Johnson RH, Norman AJ, Peacock GL, Shaver PL, Spaeth KE. 2003. Chapter 6: Livestock Nutrition, Husbandry, and Behavior. Page National Range and Pasture Handbook, 1st edition. United States Department of Agriculture Natural Resources Conservation Service. Available from http://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1043065.pdf (accessed July 11, 2016).  
> AUM_StockingRatios_byCounty.csv

Calculated area of each polygon by county, by ownership and LCC, exclusing forest, developed & urban, and strict protected areas
(overlaid PADUSCBIv2_11sgStates_nostrictPAs_byCounty_diss.shp and a masked lcc raster, based on Id field. Values in Id: first digit=LCC class; second digit=private(1) or public (0); last 5 digits=ADMIN_FIPS)
using *PerfectEnemyofGood_summariseForageAUM.r*
> /LULC/Total_area_byLCC_and_publicprivate_11states.csv

 


### Sage grouse habitat
#### PACS (priority areas for conservation) 
PACS shps were obtained from SGI & updated with data from USGS.  
PACS shp originally used has lots of holes and rogue polygons ie geometry errors USGS released a cleaned up version 2015  
> https://www.sciencebase.gov/catalog/item/560ee096e4b0ba4884c5ecf0   
> downloaded as GRSG_2015_USFWS_StatusReview_BaseData.gdb.zip  

Note the Wyoming PACS are out of date in the USGS dataset, Current (June 2016) Wyoming can be found at  
> https://wgfd.wyo.gov/WGFD/media/content/Habitat/Sage%20Grouse%20Geospatial/SageGrouseCoreAreasv4.zip  

USGS also does not have Gunnison Sage Grouse PACs.  
> downloaded as SageGrouseCoreAreasv4.zip  

*Processing notes:*  
I've combined these two datasets and updated the USGS version of the pacs with the new wyoming pacs, though not with gunnison's pacs.
Erased Wyoming from the USGS pacs (using USGS digital boundaries) then merged new wyoming pacs with usgs pacs  
> PACS_2015_USGS_USFWS_plusWyomingUpdates_Jun2016.shp   

#### Breeding habitat (lek) data 
*Source:* Doherty et al. Mapping Breeding denisities of greater sage-grouse: A tool for range-wide conservation planning.  
*Notes:* Lek locations buffered by 8.5km. Maximum count of male sage-grouse were used to map relative abundance of sage grouse. I used the breeding density_100 layer (ie. all lek locations).  

### Area of irrigated land in each county
Data on how much land is irrigated in each county was obtained from the USDA Census. This data is for harvested land (cropland) only, and excludes irrigated pasture or hay. 
*PerfectEnemyofGood_formatUSDAIrrigateddata.r*
> input: CroplandArea_Irrigation_USDACensus.csv
> output: CroplandArea_Irrigation_USDACensus_byCounty_formatted.csv
and combined with other control variables in *PerfectEnemyofGood_mergeData.r*

### Farm statistics
Statistics on farm demographics, economics etc, including areas of cropped, pasture, hay & irrigated land by county from USDA Agricultural Census 2012
> AgCensus2012.accdb
Selected statistics were formatted and combined using a query in the .accdb, 
> USDA/AgCensus/AllCountyData_AgCensus2012.csv
and combined with other control variables in *PerfectEnemyofGood_mergeData.r*

### % Urban in a county
http://www.census.gov/geo/reference/ua/urban-rural-2010.html
2010 Census Urban and Rural Classification and Urban Area Criteria

The Census Bureau’s urban-rural classification is fundamentally a delineation of geographical areas, identifying both individual urban areas and the rural areas of the nation.  The Census Bureau’s urban areas represent densely developed territory, and encompass residential, commercial, and other non-residential urban land uses.

For the 2010 Census, an urban area will comprise a densely settled core of census tracts and/or census blocks that meet minimum population density requirements, along with adjacent territory containing non-residential urban land uses as well as territory with low population density included to link outlying densely settled territory with the densely settled core.  To qualify as an urban area, the territory identified according to criteria must encompass at least 2,500 people, at least 1,500 of which reside outside institutional group quarters.  The Census Bureau identifies two types of urban areas:

Urbanized Areas (UAs) of 50,000 or more people;
Urban Clusters (UCs) of at least 2,500 and less than 50,000 people.
“Rural” encompasses all population, housing, and territory not included within an urban area.
Units seem to be m2
> USCensus_urban/PctUrbanRural_County.csv

### County population estimates
Now use annual population estimates in the models
Population estimates to 2009
https://www2.census.gov/programs-surveys/popest/tables/2000-2009/
Population estimates to 2016
https://www.census.gov/data/datasets/2016/demo/popest/counties-total.html

previously used:
Population Estimates Program   
> Downloaded from http://factfinder.census.gov/faces/nav/jsf/pages/index.xhtml
> PEP_2012_PEPTCOMP = population change  
> PEP_2012_PEPANNRES = population estimates
and combined with other control variables in *PerfectEnemyofGood_mergeData.r*  

American Community Survey  
> ACS_15_1YR_S2504 = Housing type & occupancy

### Road density
Calcuated road density for private land, by county. Roads from Jeff Evans datasets (rasters) Zone01Disturbance.gdb & Zone02Disturbance.gdb. Data is only available for western states.
Tabulate area PADUSCBIv2_Private_Land_only_AllCounty_diss.shp, then joined to this .shp
> Road_density_byCounty_privatelandsonly.csv
In R, calculated road density
Road density is proportional/ha
a$road_density=(a$roadarea/(a$roadarea+a$notroadarea))*10000
> Road_density_byCounty_privatelandsonly_CRedits.csv
and combined with other control variables in *PerfectEnemyofGood_mergeData.r*

## Combined datasets
Combined (1) Areas of private land, PACS, & sage grouse breeding habitat (2) AUMS (3) Farm statistics (4) Areas of different (NLCD) land use types (5) USDA NASS land rents (6) NRI conversion probabilities
This data was created using scripts *PerfectEnemyofGood_createInputs.r, PerfectEnemyofGood_mergeCountyData.r, JoinMosaic_SGMZ.R, JoinMosaic_SGMZ_gdbExtractor.py*
> SGCounty_alldata_160715.csv

Combined (1) JMu net returns data (NRI adjusted) (2) NASS land rents (3) Lark conversion probabilities
using *PerfectEnemyofGood_model_prelimtests.r*
> Allcounties_LarkConversion_JmuNetReturns_NASSLandrents.csv

Combined (1) Lark conversion probabilites by LCC (2) cpi adjusted NASS land rents
> Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted.csv

Combined control variables (farm demographics, popn, road density, area irrigated)
> Allcounties_ControlVariables.csv

Combined (1) Lark conversion probabilites (2) cpi adjusted NASS land rents (3) control variables
> Allcounties_LarkConversion_byLCC_NASSLandrents_CPIadjusted_plusControlVariables.csv

#Area of grassland in sgCounties (excluding mainly urban counties)
Tabulate area AnalysisAreaMask_11statesnocoast_noUrbanForest, zonal SGCounties_overlapping_PACSandBreeding_NAD83alb.shp, area is in hectares
Tabulate area AnalysisAreaMask_11statesnocoast_noUrbanForest, zonal SGCounties_private.shp, area is in hectares
> Area_grassorsagebrush_bycounty.csv
> Area_grassorsagebrush_bycounty_private.csv
In sage grouse counties, excluding CA
Total area natural vege = 1,101,653 km2
Area natural vege on private only = 428,933 km2
Median area of wet meadow mesic (1984-2011) = 50,658.3 km2


# Model notes

Model has a hard time predicting 4 points:
16077 - (Power ID) the conversion rate dramatically increased in both rate & area between 2009 & 2012 (29ha 0.0003 in 2009; 2607ha 0.034 in 2012)
30015 - (Chouteau MT) the conversion rate dramatically increased in both rate & area between 2009 & 2012 (98ha 0.0003 in 2009; 1274ha 0.004 in 2012)
56011 - (Crook WY) the conversion rate dramatically increased in both rate & area between 2009 & 2012 (434ha 0.0001 in 2009; 1492ha 0.005 in 2012)
32012 - this county (Mineral NV) has very little rangeland, consequently the conversion rate in 2012 is very large

### Zero inflation factors
tested zifs in range 1 to 1e-8. Found that ZIF of 0.0001 gave the least (flat line) skewed residuals in the Q-Q plot. 

### Mixed effects models
Trialed a series of models with range and crop rent as fixed effects, and other variables are random effects. The models had difficulty resolving when variables other than year and %crop were included. I tried splitting the data into 2 or 3 clusters, to try to estimate of betas for range rent for the different clusters of counties. These models gave significant estimates for range and crop rent, but the estimates changed dramatically depending on how clusters were chosen. In none of the models did 'cluster' fall out as a significant variable, indicating that the model was not able to significantly distinguish different betas for rent for different clusters. 

For this reason, becuase choice of clusters was so arbitrary, and influential, I chose to use a fixed effects model, where betas were estimated for all of the control variables. I performed a sensitivity analysis using a model ensemble approach where a number of models differing in control variables were run. In these set of models, range and crop rent were significant, and the estimates of betas did not differ significantly across the different models.

### Final model form
The model chosen to estimate rangeland rents was chosen from sthe set of models containing all variables,
and I performed semi-stepwise inclusion of variables, chosing the model with the lowest AICc where only significant variables are included (including year, as we know from ACR project that it is influential on rate of conversion). AICc of this model was 1409 (it doesn't make that much sense to use BIC, as we are not really interested in parsimonious models). The estimate of the coefficient for rangeland rent (beta) for this model (mod10b) is at the top end of the estimates across models.

