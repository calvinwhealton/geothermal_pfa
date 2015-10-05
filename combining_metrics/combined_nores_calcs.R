# combining risk factors except reservoirs
# potentially EGS

## reading-in the rasters
# seismic play fairway
se_3_0_3_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_3_0_3_a.tif')
se_5_0_5_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_5_0_5_a.tif')

# thermal play fairway
th_3_0_3_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_3_0_3_NA.tif')
th_5_0_5_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_5_0_5_NA.tif')

# utilization play fairway
ut5_3_0_3_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/ut5_3_0_3_NA.tif')
ut5_5_0_5_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/ut5_5_0_5_NA.tif')

## making stacked rasters
# stack for sum, product, or minimum of all
comb_pfa3_egs <- stack(c(se_3_0_3_a,th_3_0_3_NA,ut5_3_0_3_NA))
comb_pfa5_egs <- stack(c(se_5_0_5_a,th_5_0_5_NA,ut5_5_0_5_NA))

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

rm(co_5_0_5_m_egs,co_3_0_3_m_egs,co_5_0_125_p_egs,co_3_0_27_p_egs,co_3_0_9_s_egs,co_5_0_15_s_egs
   ,th_5_0_5_NA,th_3_0_3_NA,se_5_0_5_a,se_3_0_3_a,ut5_3_0_3_NA,ut5_5_0_5_NA
   ,comb_pfa3_egs,comb_pfa5_egs)
