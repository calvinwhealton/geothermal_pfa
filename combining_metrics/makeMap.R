# function to make maps for exporting to report

makeMap <- function(rast    # raster
                    ,plotnm # plot name
                    ,wd     # directory for saving file
                    ,numCol # number of colors
                    ,comTy=NA  # method of combining
                    ,numRF  # number of risk factors
                    ,sdMap=FALSE # standard deviation map T/F
                    ){
  
  # setting directory for saving file
  setwd(wd)
  
  # assigning breaks and color schemes
  if(is.na(comTy) == TRUE){
    if(numCol == 3){
      breaks <- c(-9999,-0.01,1,2,3+0.1)
      cols <- c('white','red1','yellow','green3')
    }else{
      breaks <- c(-9999,-0.01,1,2,3,4,5+0.1)
      cols <- c('white','red1','orange1','yellow','chartreuse','green3')
    }
  }else if(comTy == 'sum'){
    if(numCol == 3){
      breaks <- c(-9999,-0.01,c(1,2,3+0.1)*numRF)
      cols <- c('white','red1','yellow','green3')
    }else{
      breaks <- c(-9999,-0.01,c(1,2,3,4,5+0.1)*numRF)
      cols <- c('white','red1','orange1','yellow','chartreuse','green3')
    }
  }else if (comTy == 'prod'){
    if(numCol == 3){
      breaks <- c(-9999,-0.01,c(1,2,3+0.1)^numRF)
      cols <- c('white','red1','yellow','green3')
    }else{
      breaks <- c(-9999,-0.01,c(1,2,3,4,5+0.1)^numRF)
      cols <- c('white','red1','orange1','yellow','chartreuse','green3')
    }
  }else if(comTy == 'min'){
    if(numCol == 3){
      breaks <- c(-9999,-0.01,c(1,2,3+0.1))
      cols <- c('white','red1','yellow','green3')
    }else{
      breaks <- c(-9999,-0.01,c(1,2,3,4,5+0.1))
      cols <- c('white','red1','orange1','yellow','chartreuse','green3')
    }
  }else{
    print('Not a valid selection of comTy')
  }
  
  if(sdMap == TRUE){
    cols <- rev(brewer.pal(9,'PuBu'))
    breaks <- NA
  }
  
  # setting exporting parameters
  png(plotnm
      ,height=6
      ,width=6
      ,units='in'
      ,res=300
  )
  
  if(sdMap == TRUE){
    plot(rast
         ,col=cols # colors
         ,legend=TRUE # no legend
         ,xaxt='n' # no x axis
         ,yaxt='n' # no y axis
    )
  }else{
    # plotting raster
    plot(rast
         ,breaks=breaks
         ,col=cols # colors
         ,legend=FALSE # no legend
         ,xaxt='n' # no x axis
         ,yaxt='n' # no y axis
         )
         
    colorbar(usr)
  }
  
  # adding counties
  plot(NY_co2, add=TRUE)
  plot(PA_co2, add=TRUE)
  plot(WV_co2, add=TRUE)
  
  # adding states with thicker lines
  plot(NY2, lwd=3, add=TRUE)
  plot(PA2, lwd=3, add=TRUE)
  plot(WV2, lwd=3,add=TRUE)
    
  dev.off()
}