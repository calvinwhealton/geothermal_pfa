# function to generate histogram
# with points for thresholds

makeHist <- function(rast         # raster
                     ,thresh3     # 3-parameter threshold including min, max
                     ,thresh5     # 5-parameter threshold including min, max
                     ,rev_sc=FALSE # whether to reverse the scales (low is green)
                     ,title=''    # title of plot
                     ,xtitle=''   # label for the x-axis
                     ,xlims=c()       # limits for plotting lines
                     ,plotnm      # name of the plot
                     ,yshift=0.0025 # shift of the y location for the lines
                     ,yloc=c()        # location of the lines
                    ){
  
  # color-schemes for 3 and 5 color scheme
  # repeated values are for the lines below minimum and above maximum
  col3 <- c('red1','red1','yellow1','green3','green3')
  col5 <- c('red1','red1','orange1','yellow','chartreuse','green3','green3')

  # reversing color scheme if color scales are reversed
  if(rev_sc == TRUE){
    col3 <- rev(col3)
    col5 <- rev(col5)
  }
  
  # setting export criteria for histogram
  png(plotnm
      ,height=6
      ,width=6
      ,units='in'
      ,res=300
      )
  
  # creating histogram
  hist(rast
       ,freq=FALSE
       ,main=title
       ,xlab=xtitle
       )
  
  # setitng yloc parameter based on plot limits if not user-specified
  if(length(yloc) == 0){
    yloc <- par('usr')[3] - (par('usr')[4]-par('usr')[3])*0.25
  }

  # finding plotting limits for 3 and 5 color scheme
  if(length(xlims) == 0){
    t3 <- c(max(par("usr")[1],thresh3[1]),thresh3,min(par("usr")[2],thresh3[length(thresh3)]))
    t5 <- c(max(par("usr")[1],thresh5[1]),thresh5,min(par("usr")[2],thresh5[length(thresh5)]))
  }else{
    t3 <- c(max(par("usr")[1],xlims[1],thresh3[1]),thresh3,min(xlims[2],par("usr")[2],thresh3[length(thresh3)]))
    t5 <- c(max(par("usr")[1],xlims[1],thresh5[1]),thresh5,min(xlims[2],par("usr")[2],thresh5[length(thesh5)]))
  }
  
  
  
  # plotting lines outside of the plot region
  par(xpd=TRUE)
  
  # adding lines for the 3 color scheme
  if(length(thresh3 != 0)){
    for(i in 1:(length(thresh3)+1)){
      
      lines(c(t3[i],t3[i+1])
            ,rep(yloc-yshift,2)
            ,col=col3[i]
            ,lwd=2 
            ,lty=1)
    }
    
  }

  
  # adding lines for 5 color scheme
  if(length(thresh5) != 0){
    
    for(i in 1:(length(thresh5)+1)){
      
      lines(c(t5[i],t5[i+1])
            ,rep(yloc+yshift,2)
            ,col=col5[i]
            ,lwd=2 
            ,lty=1)
    }
    
  }

  # plotting points for 3 color scheme
  # will have outer black circle with inner white area
  if(length(thresh3 != 0)){
    points(thresh3
           ,rep(yloc-yshift,length(thresh3))
           ,col='black'
           ,pch=16
           ,cex=2)
    points(thresh3
           ,rep(yloc-yshift,length(thresh3))
           ,col='white'
           ,pch=16
           ,cex=2*0.9)
  }

  
  # plotting points for 5 color scheme
  # will have outer black circle with inner white area
  if(length(thresh5) != 0){
    
    points(thresh5
           ,rep(yloc+yshift,length(thresh5))
           ,col='black'
           ,pch=16
           ,cex=2)
    points(thresh5
           ,rep(yloc+yshift,length(thresh5))
           ,col='white'
           ,pch=16
           ,cex=2*0.9)
  }

  
  # adding text for 3 and 5 color schemes
  # text will be inside the plotted points
  if(rev_sc == FALSE){
    if(length(thresh3 != 0)){
      text(thresh3
           ,rep(yloc-yshift,length(thresh3))
           ,seq(0,3,1)
           ,cex=0.9)
    }
    if(length(thresh5 != 0)){
      text(thresh5
           ,rep(yloc+yshift,length(thresh5))
           ,seq(0,5,1)
           ,cex=0.9)
    }

  }else{
    if(length(thresh3 != 0)){
      text(thresh3
           ,rep(yloc-yshift,length(thresh3))
           ,rev(seq(0,3,1))
           ,cex=0.9)
    }

    if(length(thresh5) != 0){
      text(thresh5
           ,rep(yloc+yshift,length(thresh5))
           ,rev(seq(0,5,1))
           ,cex=0.9) 
    }

  }

  
  par(xpd=FALSE)
  dev.off()
}