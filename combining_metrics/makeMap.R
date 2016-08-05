# function to make maps for exporting to report

makeMap <- function(rast    # raster
                    ,plotnm # plot name
                    ,wd     # directory for saving file
                    ,numCol # number of colors
                    ,comTy=NA  # method of combining
                    ,numRF  # number of risk factors
                    ,sdMap=FALSE # standard deviation map T/F
                    ,leg=T # should legend be included
                    ){
  
  # setting directory for saving file
  setwd(wd)
  reclassify(rast,cbind(-9999,NA))
  temp_rast <- values(rast)
  temp_rast[which(temp_rast %in% -9999)] <- NA
  values(rast) <- temp_rast
  rast[rast < 0] <- NA
  # assigning breaks and color schemes
  if(is.na(comTy) == TRUE){
    if(numCol == 3){
      breaks <- c(0,1,2,3)
      cols <- c('red1','yellow','forestgreen')
    }else{
      breaks <- c(0,1,2,3,4,5)
      cols <- c('red1','orange1','yellow','chartreuse','forestgreen')
    }
  }else if(comTy == 'sum'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3)*numRF)
      cols <- c('red1','yellow','forestgreen')
    }else{
      breaks <- c(0,c(1,2,3,4,5+0.1)*numRF)
      cols <- c('red1','orange1','yellow','chartreuse','forestgreen')
    }
  }else if (comTy == 'prod'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3)^numRF)
      cols <- c('red1','yellow','forestgreen')
    }else{
      breaks <- c(0,c(1,2,3,4,5)^numRF)
      cols <- c('red1','orange1','yellow','chartreuse','forestgreen')
    }
  }else if(comTy == 'min'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3))
      cols <- c('red1','yellow','forestgreen')
    }else{
      breaks <- c(0,c(1,2,3,4,5))
      cols <- c('red1','orange1','yellow','chartreuse','forestgreen')
    }
  }else{
    print('Not a valid selection of comTy')
  }
  
  if(sdMap == TRUE){
    cols <- rev(brewer.pal(9,'GnBu')[seq(3,9)])
    breaks <- round(seq(0,max(values(rast),na.rm=T),length.out=length(cols)+1),digits=2)
  }
  
  # setting exporting parameters
#   pdf('test.pdf'
#       ,height=6
#       ,width=6)
  png(plotnm
      ,height=6
      ,width=6
      ,units='in'
      ,res=300
  )
  
  if(sdMap == TRUE){
    raster:plot(rast
                ,breaks=breaks
                ,col=cols # colors
                #,legend=leg # include legend
                #,xaxt='n' # no x axis
                #,yaxt='n' # no y axis
                ,ext=extent(380000,1050000,4100000,4850000)
                #,ylim=c(4150000,4850000)
                #,xlim=c(400000,1000000)
                #,axes=F
            
)
  }else{
    
    # plotting raster
    raster:plot(rast
         ,breaks=breaks
         ,col=cols # colors
         #,legend=leg # include legend
         #,xaxt='n' # no x axis
         #,yaxt='n' # no y axis
         ,ext=extent(380000,1050000,4100000,4850000)
         #,ylim=c(4150000,4850000)
         #,xlim=c(400000,1000000)
         #,axes=F
         )
  }
  #axis(2,at=c(41,43,45,47)*10^5,labels=c(41,43,45,47)*10^5)
  #axis(1,at=c(4,6,8,10)*10^5,labels=format(c(4,6,8,10)*10^5,scientific=F))
  
  #map.axes()
       
       
  # adding counties
  plot(NY_co2, add=TRUE)
  plot(PA_co2, add=TRUE)
  plot(WV_co2, add=TRUE)
  
  # adding states with thicker lines
  plot(NY2, lwd=3, add=TRUE)
  plot(PA2, lwd=3, add=TRUE)
  plot(WV2, lwd=3,add=TRUE)
 
  addnortharrow(pos='topleft',padin=c(0.15,0.15),scale=0.5)
  addscalebar(plotunit='m',padin=c(0.15,0.7),pos='topleft')
  
  text(x=800000,y=4200000,'Projection: UTM 17N')
  dev.off()
}