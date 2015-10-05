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

##### defining working directories #####
wd_raster <- '/Users/calvinwhealton/Dropbox/PFA_Rasters'
wd_image <- '/Users/calvinwhealton/Desktop/combining_rf'
wd_source <- '/Users/calvinwhealton/GitHub/geothermal/combining_metrics'

##### importing functions/libraries #####
setwd(wd_source)
source('libraries_funcs.R')

##### loading-in state/county shapefiles #####
setwd(wd_source)
source('shapefiles_rasters_general.R')

##### THERMAL ######
# file containing all commands for analysis of thermal resource
setwd(wd_source)
source('therm_d80_calcs.R') 

##### RESERVOIR ######
setwd(wd_source)
source('reservoir_calcs.R')

##### UTILIZATION ####
setwd(wd_source)
source('utilization_calcs.R')

##### SEISMIC ######
setwd(wd_source)
source('seismic_calcs.R')

#### combining the metrics ####
# for all risk factors
setwd(wd_source)
source('combine_all_calcs.R')

# for geology risk factors
setwd(wd_source)
source('combined_geology_calcs.R')

# for no reservoirs (EGS) risk factors
setwd(wd_source)
source('combined_nores_calcs.R')


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




