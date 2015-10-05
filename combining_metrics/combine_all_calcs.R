# script for combining all rasters and
# creating output plots

## reading-in the rasters
# seismic play fairway
se_3_0_3_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_3_0_3_a.tif')
se_5_0_5_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_5_0_5_a.tif')

a <- values(se_3_0_3_a)
a[(a == -9999)] <- NA
values(se_3_0_3_a) <- a

se_3_0_3_a@file@nodatavalue <- NA

# reservoir play fairway
re_3_0_3_x <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_3_0_3_x.tif')
re_5_0_5_x <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_5_0_5_x.tif')

# thermal play fairway
th_3_0_3_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_3_0_3_NA.tif')
th_5_0_5_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_5_0_5_NA.tif')

# utilization play fairway
ut5_3_0_3_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/ut5_3_0_3_NA.tif')
ut5_5_0_5_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/ut5_5_0_5_NA.tif')

## making stacked rasters
# stack for sum, product, or minimum of all
comb_pfa3_all <- stack(c(se_3_0_3_a,re_3_0_3_x,th_3_0_3_NA,ut5_3_0_3_NA))
comb_pfa5_all <- stack(c(se_5_0_5_a,re_5_0_5_x,th_5_0_5_NA,ut5_5_0_5_NA))


## calculating the combined metrics
# all variables
co_3_0_12_s <- calc(comb_pfa3_all,fun=sum,na.rm=FALSE)
co_5_0_20_s <- calc(comb_pfa5_all,fun=sum,na.rm=FALSE)

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
co_3_0_81_p <- calc(comb_pfa3_all,fun=prod,na.rm=FALSE)
co_5_0_625_p <- calc(comb_pfa5_all,fun=prod,na.rm=FALSE)

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
co_3_0_3_m <- calc(comb_pfa3_all,fun=min,na.rm=FALSE)
co_5_0_5_m <- calc(comb_pfa5_all,fun=min,na.rm=FALSE)

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

# deleting files
rm(co_5_0_5_m,co_3_0_3_m,co_5_0_20_s,co_3_0_12_s,co_5_0_625_p,co_3_0_81_p
   ,ut5_5_0_5_NA,ut5_3_0_3_NA,th_5_0_5_NA,th_3_0_3_NA
   ,se_5_0_5_a,se_3_0_3_a,re_5_0_5_x,re_3_0_3_x
   ,comb_pfa3_all,comb_pfa5_all)