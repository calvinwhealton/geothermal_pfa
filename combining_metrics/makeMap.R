# function to make maps for exporting to report

makeMap <- function(rast    # raster
                    ,plotnm # plot name
                    ,wd     # directory for saving file
                    ,numCol # number of colors
                    ,comTy=NA  # method of combining
                    ,numRF  # number of risk factors
                    ,sdMap=FALSE # standard deviation map T/F
                    ,leg2=F # should legend include thresholds only (FALSE), or also unscaled values (TRUE)
                    ,grey=F #Grey scale maps
                    ,County=F #include counties
                    ,RawThreshVals=NA #Vector of unscaled values corresponding to the thresholds. Only used to plot these on the legend.
                    ,Unit = NA #Unit for the RawThreshVals
                    ,dpi = 300 #dots per square inch for the figure
                    ,FigFun = 'png' #figure function to use. Accepts tiff and png
                    ){
  
  # setting directory for saving file
  setwd(wd)
  reclassify(rast,cbind(-9999,NA))
  temp_rast <- values(rast)
  temp_rast[which(temp_rast %in% -9999)] <- NA
  values(rast) <- temp_rast
  rast[rast < 0] <- NA
  
  # assigning color schemes
  if (grey){
    cols <- rev(grey.colors(numCol))
  }else{
    if(numCol == 3){
      cols <- c('red1','yellow','forestgreen')
    }
    else{
      cols <- c('red1','orange1','yellow','chartreuse','forestgreen')
    }
  }
  
  # assigning breaks 
  if(is.na(comTy) == TRUE){
    if(numCol == 3){
      breaks <- c(0,1,2,3)
    }else{
      breaks <- c(0,1,2,3,4,5)
    }
  }else if(comTy == 'sum'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3)*numRF)
    }else{
      breaks <- c(0,c(1,2,3,4,5+0.1)*numRF)
    }
  }else if (comTy == 'prod'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3)^numRF)
    }else{
      breaks <- c(0,c(1,2,3,4,5)^numRF)
    }
  }else if(comTy == 'min'){
    if(numCol == 3){
      breaks <- c(0,c(1,2,3))
    }else{
      breaks <- c(0,c(1,2,3,4,5))
    }
  }else{
    print('Not a valid selection of comTy')
  }
  
  if(sdMap == TRUE){
    if (grey){
      breaks <- pretty(seq(0,max(values(rast), na.rm=T), length.out = 6))
      cols <- rev(grey.colors(length(breaks)))
    }else{
      cols <- rev(brewer.pal(9,'GnBu')[seq(3,9)])
      breaks <- pretty(seq(0,max(values(rast),na.rm=T),length.out=length(cols)+1))
    }
  }
  
  # setting exporting parameters
#   pdf('test.pdf'
#       ,height=6
#       ,width=6)
  if (FigFun == 'png'){
    png(plotnm
        ,height=6
        ,width=6
        ,units='in'
        ,res=dpi
    )
  }else{
    tiff(plotnm
        ,height=6
        ,width=6
        ,units='in'
        ,res=dpi
    )
  }
  
  
  if(sdMap == TRUE){
    plot(rast
         ,breaks=breaks
         ,col=cols # colors
        #,legend=leg # include legend
        #,xaxt='n' # no x axis
        #,yaxt='n' # no y axis
         ,ext=extent(350000,1050000,4100000,4850000)
        #,ylim=c(4150000,4850000)
        #,xlim=c(400000,1000000)
        ,axes=F
        )
    
    if (County){
      # adding counties
      plot(NY_co2, add=TRUE)
      plot(PA_co2, add=TRUE)
      plot(WV_co2, add=TRUE)
    }     
    
    # adding states with thicker lines
    plot(NY2, lwd=3, add=TRUE)
    plot(PA2, lwd=3, add=TRUE)
    plot(WV2, lwd=3,add=TRUE)
    
    addnortharrow(pos='topleft',padin=c(0.15,0.15),scale=0.5)
    addscalebar(plotunit='m',padin=c(0.15,0.7),pos='topleft')
    
    text(x=800000,y=4200000,'Projection: UTM 17N')
    axis(2,at=c(42,44,46,48)*10^5,labels=format(c(42,44,46,48)*10^5,scientific=T))
    axis(1,at=c(4,6,8,10)*10^5,labels=format(c(4,6,8,10)*10^5,scientific=T))
  }else{
    if (leg2){
      # plotting raster with 2 legends

      #First plot map with no legend
      plot(rast
           ,breaks=breaks
           ,col=cols # colors
           ,legend = FALSE
           #,xaxt='n' # no x axis
           #,yaxt='n' # no y axis
           ,ext=extent(350000,1050000,4100000,4850000)
           #,ylim=c(4150000,4850000)
           #,xlim=c(400000,1000000)
           ,axes=F
      )
      
      if (County){
        # adding counties
        plot(NY_co2, add=TRUE)
        plot(PA_co2, add=TRUE)
        plot(WV_co2, add=TRUE)
      }     
      
      # adding states with thicker lines
      plot(NY2, lwd=3, add=TRUE)
      plot(PA2, lwd=3, add=TRUE)
      plot(WV2, lwd=3,add=TRUE)
      
      addnortharrow(pos='topleft',padin=c(0.15,0.15),scale=0.5)
      addscalebar(plotunit='m',padin=c(0.15,0.7),pos='topleft')
      
      text(x=800000,y=4200000,'Projection: UTM 17N')
      axis(2,at=c(42,44,46,48)*10^5,labels=format(c(42,44,46,48)*10^5,scientific=T))
      axis(1,at=c(4,6,8,10)*10^5,labels=format(c(4,6,8,10)*10^5,scientific=T))
      
      #Add Unit to plot:
      mtext(text = Unit, side = 1, at = 12e5, cex = 1.1, line = -18)

      #Then plot the threshold legend without the map.
      plot(rast
           ,breaks=breaks
           ,col=cols # colors
           ,legend.only = TRUE
           #,xaxt='n' # no x axis
           #,yaxt='n' # no y axis
           ,ext=extent(350000,1050000,4100000,4850000)
           #,ylim=c(4150000,4850000)
           #,xlim=c(400000,1000000)
           ,axes=F
           ,axis.args = list(at = breaks, labels = breaks, mgp = c(-3,-1.5,0), tck = 1.7)
      )

      #Then plot the unscaled legend without the map.
      plot(rast
           ,breaks=breaks
           ,col=cols # colors
           ,legend.only = TRUE
           #,xaxt='n' # no x axis
           #,yaxt='n' # no y axis
           ,ext=extent(350000,1050000,4100000,4850000)
           #,ylim=c(4150000,4850000)
           #,xlim=c(400000,1000000)
           ,axes=F
           ,axis.args = list(at = seq(0,5,1), labels = RawThreshVals, mgp = c(3,1,0))
      )
    }else{
      # plotting raster
      plot(rast
           ,breaks=breaks
           ,col=cols # colors
           #,legend=leg # include legend
           #,xaxt='n' # no x axis
           #,yaxt='n' # no y axis
           ,ext=extent(350000,1050000,4100000,4850000)
           #,ylim=c(4150000,4850000)
           #,xlim=c(400000,1000000)
           ,axes=F
      )
      
      if (County){
        # adding counties
        plot(NY_co2, add=TRUE)
        plot(PA_co2, add=TRUE)
        plot(WV_co2, add=TRUE)
      }     
      
      # adding states with thicker lines
      plot(NY2, lwd=3, add=TRUE)
      plot(PA2, lwd=3, add=TRUE)
      plot(WV2, lwd=3,add=TRUE)
      
      addnortharrow(pos='topleft',padin=c(0.15,0.15),scale=0.5)
      addscalebar(plotunit='m',padin=c(0.15,0.7),pos='topleft')
      
      text(x=800000,y=4200000,'Projection: UTM 17N')
      axis(2,at=c(42,44,46,48)*10^5,labels=format(c(42,44,46,48)*10^5,scientific=T))
      axis(1,at=c(4,6,8,10)*10^5,labels=format(c(4,6,8,10)*10^5,scientific=T))
    }
  }
  
  #map.axes()
  dev.off()
}