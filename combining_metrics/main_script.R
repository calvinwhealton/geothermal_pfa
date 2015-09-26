##################################################
# script for combining risk factor rasters (RFs)
# into the play fairway metric raster
#
# steps include:
# -1. importing functions/libraries
# 0. importing rasters
# 1. verifying that the rasters are in the same system
# 2. converting individual rasters into the play fairway rank
# 3. combining the RFs using a specified function
# 4. saving output

##### importing functions/libraries #####

# libraries
library(sp)
library(raster)
library(rgdal)
library(rasterVis)
library(maps)
library(maptools)
library(xlsx)
library(pracma)
library(RColorBrewer)

##### defining working directories #####
wd_raster <- '/Users/calvinwhealton/Dropbox/PFA_Rasters'
wd_image <- '/Users/calvinwhealton/Desktop/combining_rf'
wd_source <- '/Users/calvinwhealton/GitHub/geothermal/combining_metrics'

##### loading-in state/county shapefiles #####
# states
States = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# counties
Counties = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_county_WGS84_prj", stringsAsFactors=FALSE)
NY_co = Counties[which(Counties$STATEFP == "36"),]
PA_co = Counties[which(Counties$STATEFP == "42"),]
WV_co = Counties[which(Counties$STATEFP == "54"),]

# converting to UTM 17N system
NY2 <- spTransform(NY,CRS("+init=epsg:31986"))
PA2 <- spTransform(PA,CRS("+init=epsg:31986"))
WV2 <- spTransform(WV,CRS("+init=epsg:31986"))

NY_co2 <- spTransform(NY_co,CRS("+init=epsg:31986"))
PA_co2 <- spTransform(PA_co,CRS("+init=epsg:31986"))
WV_co2 <- spTransform(WV_co,CRS("+init=epsg:31986"))

# importing city locations
cities <- readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="usCensusPlaces", stringsAsFactors=FALSE)
cities2 <- spTransform(cities,CRS("+init=epsg:31986"))

# importing grid locations
fishnet <- readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="Fishnet2_label", stringsAsFactors=FALSE)



##### user-defined functions #####
setwd('/Users/calvinwhealton/GitHub/geothermal/combining_metrics')
source('checkSameProjCoords.R')
source('convertRasterPFAMetric.R')
source('combineRFs.R')
source('makeUtilBuf.R')
source('makeHist.R')
source('makeMap.R')
source('saveRast.R')
source('plotWeightBuf.R')

##### THERMAL ######
# file containing all commands for analysis of thermal resource
setwd(wd_source)
source('thermal_calcs.R') 



##### RESERVOIR ######
setwd(wd_source)
source('reservoir_calcs.R')


##### UTILIZATION ####
setwd(wd_source)
source('utilization_calcs.R')

##### SEISMIC ######
setwd(wd_source)
source('seismic_calcs.R')

# using sums
co_3_0_12_s <- calc(comb_pfa3,fun=sum,na.rm=FALSE)
co_5_0_20_s <- calc(comb_pfa5,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_12_s
         ,thresh3=c(0,1,2,3)*4
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_12_s_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_12_s
         ,wd=wd_raster
         ,rastnm='co_3_0_12_s.tif')
makeMap(rast=co_3_0_12_s
        ,plotnm='co_3_0_12_s.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_20_s
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*4
         ,rev_sc=FALSE
         ,plotnm='co_5_0_20_s_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_20_s
         ,wd=wd_raster
         ,rastnm='co_5_0_20_s.tif')
makeMap(rast=co_5_0_20_s
        ,plotnm='co_5_0_20_s.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=4)

# using products
co_3_0_81_p <- calc(comb_pfa3,fun=prod,na.rm=FALSE)
co_5_0_625_p <- calc(comb_pfa5,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_81_p
         ,thresh3=c(0,1,2,3)^4
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_81_p_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_81_p
         ,wd=wd_raster
         ,rastnm='co_3_0_81_p.tif')
makeMap(rast=co_3_0_81_p
        ,plotnm='co_3_0_81_p.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_625_p
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^4
         ,rev_sc=FALSE
         ,plotnm='co_5_0_625_p_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_625_p
         ,wd=wd_raster
         ,rastnm='co_5_0_625_p.tif')
makeMap(rast=co_5_0_625_p
        ,plotnm='co_5_0_625_p.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=4)

# using minimums
co_3_0_3_m <- calc(comb_pfa3,fun=min,na.rm=FALSE)
co_5_0_5_m <- calc(comb_pfa5,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_hist.png'
         ,yloc=-1.2
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_3_m
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m.tif')
makeMap(rast=co_3_0_3_m
        ,plotnm='co_3_0_3_m.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_5_m
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_5_m
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m.tif')
makeMap(rast=co_5_0_5_m
        ,plotnm='co_5_0_5_m.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=1)


# geologic only
comb_pfa3_geo <- stack(c(re_3_0_3_NA,th_3_0_3_NA,se_3_0_3_a))
comb_pfa5_geo <- stack(c(re_5_0_5_NA,th_5_0_5_NA,se_5_0_5_a))

# sums
co_3_0_9_s_geo <- calc(comb_pfa3_geo,fun=sum,na.rm=FALSE)
co_5_0_15_s_geo <- calc(comb_pfa5_geo,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_9_s_geo
         ,thresh3=c(0,1,2,3)*3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_9_s_geo_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_9_s_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_9_s_geo.tif')
makeMap(rast=co_3_0_9_s_geo
        ,plotnm='co_3_0_9_s_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_15_s_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_15_s_geo_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_15_s_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_15_s_geo.tif')
makeMap(rast=co_5_0_15_s_geo
        ,plotnm='co_5_0_15_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=3)

# product
co_3_0_27_p_geo <- calc(comb_pfa3_geo,fun=prod,na.rm=FALSE)
co_5_0_125_p_geo <- calc(comb_pfa5_geo,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_27_p_geo
         ,thresh3=c(0,1,2,3)^3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_27_p_geo_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_27_p_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_27_p_geo.tif')
makeMap(rast=co_3_0_27_p_geo
        ,plotnm='co_3_0_27_p_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_125_p_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_125_p_geo_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_125_p_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_125_p_geo.tif')
makeMap(rast=co_5_0_125_p_geo
        ,plotnm='co_5_0_125_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=3)

# minimum
co_3_0_3_m_geo <- calc(comb_pfa3_geo,fun=min,na.rm=FALSE)
co_5_0_5_m_geo <- calc(comb_pfa5_geo,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m_geo
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_geo_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_3_m_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m_geo.tif')
makeMap(rast=co_3_0_3_m_geo
        ,plotnm='co_3_0_3_m_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_5_m_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_geo_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_5_m_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m_geo.tif')
makeMap(rast=co_5_0_5_m_geo
        ,plotnm='co_5_0_5_m_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3)

# no reservoirs
comb_pfa3_egs <- stack(c(ut5_3_0_3_NA,th_3_0_3_NA,se_3_0_3_a))
comb_pfa5_egs <- stack(c(ut5_5_0_5_NA,th_5_0_5_NA,se_5_0_5_a))

# sums
co_3_0_9_s_egs <- calc(comb_pfa3_egs,fun=sum,na.rm=FALSE)
co_5_0_15_s_egs <- calc(comb_pfa5_egs,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_9_s_egs
         ,thresh3=c(0,1,2,3)*3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_9_s_egs_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_9_s_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_9_s_egs.tif')
makeMap(rast=co_3_0_9_s_egs
        ,plotnm='co_3_0_9_s_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_15_s_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_15_s_egs_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_15_s_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_15_s_egs.tif')
makeMap(rast=co_5_0_15_s_egs
        ,plotnm='co_5_0_15_s_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=3)

# product
co_3_0_27_p_egs <- calc(comb_pfa3_egs,fun=prod,na.rm=FALSE)
co_5_0_125_p_egs <- calc(comb_pfa5_egs,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_27_p_egs
         ,thresh3=c(0,1,2,3)^3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_27_p_egs_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_27_p_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_27_p_egs.tif')
makeMap(rast=co_3_0_27_p_egs
        ,plotnm='co_3_0_27_p_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_125_p_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_125_p_egs_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_125_p_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_125_p_egs.tif')
makeMap(rast=co_5_0_125_p_egs
        ,plotnm='co_5_0_125_p_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=3)

# minimum
co_3_0_3_m_egs <- calc(comb_pfa3_egs,fun=min,na.rm=FALSE)
co_5_0_5_m_egs <- calc(comb_pfa5_egs,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m_egs
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_egs_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_3_m_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m_egs.tif')
makeMap(rast=co_3_0_3_m_egs
        ,plotnm='co_3_0_3_m_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_5_m_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_egs_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_5_m_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m_egs.tif')
makeMap(rast=co_5_0_5_m_egs
        ,plotnm='co_5_0_5_m_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3)


########################
# extracting values of layers for cities
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA,se_5_0_5_a,ut5_5_0_5_NA
                                ,co_5_0_20_s,co_5_0_625_p,co_5_0_5_m))
                
combined_all_rasters[combined_all_rasters<0] <- NA
extract_pts <- extract(x=combined_all_rasters
              ,y=cities2
              ,sp=TRUE
              ,nl=7
              ,df=TRUE
              ,na.rm=TRUE
              ,method='simple'
              ,buffer=2000
              ,fun=mean
              )

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[15] <- "th_5"
names(extract_pts2)[16] <- "re_5"
names(extract_pts2)[17] <- "se_5"
names(extract_pts2)[18] <- "ut_5"
names(extract_pts2)[19] <- "co_s_5"
names(extract_pts2)[20] <- "co_p_5"
names(extract_pts2)[21] <- "co_m_5"


extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]

ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5 > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5 > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
png('scatter_p_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]
     ,extract_pts3$co_p_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Sum'
     ,ylab='Combined Product'
     )
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Sum'
     ,ylab='Combined Minimum'
)
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_p.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_p_5[inds_gtr0]
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Product'
     ,ylab='Combined Minimum'
)
dev.off()

inds_top10pct <- which(extract_pts3$co_s_5 > quantile(extract_pts3$co_s_5,0.9))

# making plot with color
setwd(wd_image)
png('parallel.png'
    ,height=4
    ,width=6
    ,units='in'
    ,res=300
)

plot(NA,NA
     ,xlim=c(0,3)
     ,ylim=c(0,5)
     ,xaxt='n'
     ,xlab=''
     ,ylab='Scaled Risk Factor')
lines(c(1,1)
     ,c(-1,6)
     ,lwd=1
     ,col='black')
lines(c(2,2)
      ,c(-1,6)
      ,lwd=1
      ,col='black')
lines(c(3,3)
      ,c(-1,6)
      ,lwd=1
      ,col='black')

for(i in 1:length(inds_top10pct)){
  
  lines(c(0,1,2,3)
       ,c(extract_pts3$th_5[inds_top10pct[i]],extract_pts3$re_5[inds_top10pct[i]],extract_pts3$se_5[inds_top10pct[i]],extract_pts3$ut_5[inds_top10pct[i]])
       ,lwd=2
       ,col='seagreen')
  
}
par(xpd=TRUE)
text(c(0,1,2,3)
     ,y=-0.5
     ,labels=c('Thermal','Reservoir', 'Seismic','Utilization')
     ,col='black'
     ,adj=0.5)
par(xpd=FALSE)
dev.off()




