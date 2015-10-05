# script for combining geologic risk factors only
# and approximate uncertainty

## reading-in the rasters
# seismic play fairway
se_3_0_3_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_3_0_3_a.tif')
se_5_0_5_a <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_5_0_5_a.tif')

# seismic uncertainty
se_pfa_var3 <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_pfa_var3.tif')
se_pfa_var5 <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_pfa_var5.tif')

# reservoir play fairway
re_3_0_3_m <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_3_0_3_x.tif')
re_5_0_5_x <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_5_0_5_x.tif')

# reservoir uncertainty
re_3_0_3_x <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_pfa_var3.tif')
re_5_0_5_x <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_pfa_var5.tif')

# thermal play fairway
th_3_0_3_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_3_0_3_NA.tif')
th_5_0_5_NA <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_5_0_5_NA.tif')

# thermal uncertainty
th_pfa_var3 <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_pfa_var3.tif')
th_pfa_var5 <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_pfa_var5.tif')

## making stacked rasters
# stack for sum, product, or minimum of geology (no utilization)
comb_pfa3_geo <- stack(c(se_3_0_3_a,re_3_0_3_x,th_3_0_3_NA))
comb_pfa5_geo <- stack(c(se_5_0_5_a,re_5_0_5_x,th_5_0_5_NA))

comb_pfa3_geo_sum_uncer <- stack(c(se_pfa_var3,re_pfa_var3,th_pfa_var3))
comb_pfa5_geo_sum_uncer <- stack(c(se_pfa_var5,re_pfa_var5,th_pfa_var5))

# stacking two of the last elements because the values need to be squared for individual terms
comb_pfa3_geo_prod_uncer1 <- stack(c(se_pfa_var3,re_3_0_3_x,re_3_0_3_x,th_3_0_3_NA,th_3_0_3_NA))
comb_pfa3_geo_prod_uncer2 <- stack(c(re_pfa_var3,se_3_0_3_a,se_3_0_3_a,th_3_0_3_NA,th_3_0_3_NA))
comb_pfa3_geo_prod_uncer3 <- stack(c(th_pfa_var3,se_3_0_3_a,se_3_0_3_a,re_3_0_3_x,re_3_0_3_x))

comb_pfa5_geo_prod_uncer1 <- stack(c(se_pfa_var5,re_5_0_5_x,re_5_0_5_x,th_5_0_5_NA,th_5_0_5_NA))
comb_pfa5_geo_prod_uncer2 <- stack(c(re_pfa_var5,se_5_0_5_a,se_5_0_5_a,th_5_0_5_NA,th_5_0_5_NA))
comb_pfa5_geo_prod_uncer3 <- stack(c(th_pfa_var5,se_5_0_5_a,se_5_0_5_a,re_5_0_5_x,re_5_0_5_x))

# calcuting the certainty for each term in the Taylor series expansion
# var(a*b*c) ~= var(a)*(bc)^2 + var(b)*(ac)^2 + var(c)*(ab)^2
comb_pfa3_geo_prod_t1 <- calc(comb_pfa3_geo_prod_uncer1,fun=prod,rm.na=FALSE)
comb_pfa3_geo_prod_t2 <- calc(comb_pfa3_geo_prod_uncer2,fun=prod,rm.na=FALSE)
comb_pfa3_geo_prod_t3 <- calc(comb_pfa3_geo_prod_uncer3,fun=prod,rm.na=FALSE)

comb_pfa5_geo_prod_t1 <- calc(comb_pfa5_geo_prod_uncer1,fun=prod,rm.na=FALSE)
comb_pfa5_geo_prod_t2 <- calc(comb_pfa5_geo_prod_uncer2,fun=prod,rm.na=FALSE)
comb_pfa5_geo_prod_t3 <- calc(comb_pfa5_geo_prod_uncer3,fun=prod,rm.na=FALSE)

# geologic only
comb_pfa3_geo <- stack(c(re_3_0_3_x,th_3_0_3_NA,se_3_0_3_a))
comb_pfa5_geo <- stack(c(re_5_0_5_x,th_5_0_5_NA,se_5_0_5_a))

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

## uncertainty calculations
# product
comb_pfa3_geo_prod_var <- calc(stack(c(comb_pfa3_geo_prod_t1,comb_pfa3_geo_prod_t2,comb_pfa3_geo_prod_t3)),fun=sum,na.rm=FALSE)
comb_pfa5_geo_prod_var <- calc(stack(c(comb_pfa5_geo_prod_t1,comb_pfa5_geo_prod_t2,comb_pfa5_geo_prod_t3)),fun=sum,na.rm=FALSE)

comb_pfa3_geo_sum_var <- calc(stack(c(re_pfa_var3,th_pfa_var3,se_pfa_var3)),fun=sum,na.rm=FALSE)
comb_pfa5_geo_sum_var <- calc(stack(c(re_pfa_var5,th_pfa_var5,se_pfa_var5)),fun=sum,na.rm=FALSE)


saveRast(rast=comb_pfa3_geo_prod_var
         ,wd=wd_raster
         ,rastnm='co_var_3_geo_prod.tif')
makeMap(rast=comb_pfa3_geo_prod_var
        ,plotnm='co_var_3_geo_prod.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=comb_pfa5_geo_prod_var
         ,wd=wd_raster
         ,rastnm='co_var_5_geo_prod.tif')
makeMap(rast=comb_pfa5_geo_prod_var
        ,plotnm='co_var_5_geo_prod.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=comb_pfa3_geo_sum_var
         ,wd=wd_raster
         ,rastnm='co_var_3_geo_sum.tif')
makeMap(rast=comb_pfa3_geo_sum_var
        ,plotnm='co_var_3_geo_sum.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=comb_pfa3_geo_prod_var
         ,wd=wd_raster
         ,rastnm='co_var_5_geo_sum.tif')
makeMap(rast=comb_pfa3_geo_prod_var
        ,plotnm='co_var_5_geo_sum.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)



saveRast(rast=calc(comb_pfa3_geo_prod_var,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_sd_3_geo_prod.tif')
makeMap(rast=calc(comb_pfa3_geo_prod_var,fun=sqrt)
        ,plotnm='co_sd_3_geo_prod.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(comb_pfa5_geo_prod_var,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_sd_5_geo_prod.tif')
makeMap(rast=calc(comb_pfa5_geo_prod_var,fun=sqrt)
        ,plotnm='co_sd_5_geo_prod.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(comb_pfa3_geo_sum_var,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_sd_3_geo_sum.tif')
makeMap(rast=calc(comb_pfa3_geo_sum_var,fun=sqrt)
        ,plotnm='co_sd_3_geo_sum.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(comb_pfa5_geo_sum_var,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_sd_5_geo_sum.tif')
makeMap(rast=calc(comb_pfa5_geo_sum_var,fun=sqrt)
        ,plotnm='co_sd_5_geo_sum.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

rm(re_3_0_3_x,re_5_0_5_x,se_3_0_3_a,se_5_0_5_a,th_3_0_3_NA,th_5_0_5_NA,
   comb_pfa5_geo_sum_var,comb_pfa3_geo_sum_var,comb_pfa5_geo_prod_var,comb_pfa3_geo_prod_var,
   re_pfa_var3,re_pfa_var5,th_pfa_var3,th_pfa_var5,se_pfa_var3,se_pfa_var5,
   comb_pfa3_geo_sum_uncer,comb_pfa5_geo_sum_uncer,
   comb_pfa3_prod_uncer1,comb_pfa3_prod_uncer2,comb_pfa3_prod_uncer3,
   comb_pfa3_prod_t1,comb_pfa3_prod_t2,comb_pfa3_prod_t3)

