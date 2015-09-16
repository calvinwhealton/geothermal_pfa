# function to make image of weighting matrix

makeWeightBuf <- function(dist # distance of buffer
                          ,wd # working directory
                          ,plotnm # plot name
                          ){
  
  # setting working directory
  setwd(wd)
  
  # setting exporting parameters for plot
  png(plotnm
      ,height=6
      ,width=6
      ,units='in'
      ,res=300
  )
  
  # making buffer
  image(makeUtilBufWeight(dist) # calculating weighting matrix
        ,xaxt='n' # no x axis labels
        ,yaxt='n' # no y axis labels
  )
  
  # adding vertical lies
  for(i in 0:(2*dist)){
    lines(c(i+0.5,i+0.5)/(2*dist),c(-0.5,1.5)
          ,lty=2 # line type
          ,lwd=2 # line weight
          ,col='gray48' # color
          )
  }
  
  # adding horizontal lines
  for(i in 0:(2*dist)){
    lines(c(-0.5,1.5),c(i+0.5,i+0.5)/(2*dist)
          ,lty=2 # line type
          ,lwd=2 # line weight
          ,col='gray48' # color
          )
  }
  
  # adding center location
  points(0.5
         ,0.5
        ,pch=19
        ,col='black'
        ,cex=1.2)
  
  
  dev.off()
}