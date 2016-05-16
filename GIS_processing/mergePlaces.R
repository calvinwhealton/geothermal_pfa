###################
# script for merging US Census Places
# main steps are:
# 1. loading libraries
# 2. loading US Census places shapefile
# 3. determining which places are candidates for merging
# 4. mering values
# 5. repeating 3 and 4 until no more merges are needed
# 6. saving merged shapefile

# loading libraries
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
library(spdep)        # other spatial processing
library(rgeos)        # for buffering places of interest and finding intersections
library(maptools)
library()
#library(xlsx)         # for reading in xlsx tables
#library(RColorBrewer) # R color brewer palettes for parallel axis plot
#library(pracma)       # for interpolation in tables
#library(vioplot)      # for violin plots


# reading-in shapefile with US Census "places" clipped to states and region
setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/cb_2015_us_state_500k_copy_clipAppSt')
places <- readOGR(dsn='/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/cb_2015_us_state_500k_copy_clipAppSt'
                 ,layer='clipped_appBas_st')

# transforming to UTM 17N
# units in this projection are in meters
placesTrans <- spTransform(places,CRS("+init=epsg:31986"))

# buffering the shapefile polygons by a small amount
# accounts for polygons might not intersect exactly
# so gTouches() would not work quite as well
placesTransTouch <- gTouches(placesTrans)

# buffering polygons by 25 m
placesTransBuf25 <- gBuffer(placesTrans # shapefile
                            ,byid=T
                            ,width=25 # in units of projection
)

# buffering polygons by 50 m
placesTransBuf50 <- gBuffer(placesTrans # shapefile
                             ,byid=T
                             ,width=50 # in units of projection
)

# buffering polygons by 75 m
placesTransBuf75 <- gBuffer(placesTrans # shapefile
                            ,byid=T
                            ,width=75 # in units of projection
)

# buffering polygons by 100 m
placesTransBuf100 <- gBuffer(placesTrans # shapefile
                              ,byid=T
                              ,width=100 # in units of projection
                              )

# buffering polygons by 125 m
placesTransBuf125 <- gBuffer(placesTrans # shapefile
                             ,byid=T
                             ,width=125 # in units of projection
)

# buffering polygons by 150 m
placesTransBuf150 <- gBuffer(placesTrans # shapefile
                             ,byid=T
                             ,width=150 # in units of projection
)

# buffering polygons by 175 m
placesTransBuf175 <- gBuffer(placesTrans # shapefile
                             ,byid=T
                             ,width=175 # in units of projection
)

# buffering polygons by 200 m
placesTransBuf200 <- gBuffer(placesTrans # shapefile
                             ,byid=T
                             ,width=200 # in units of projection
)


placesTransCont25 <- gIntersects(placesTransBuf25,byid = T)
placesTransCont50 <- gIntersects(placesTransBuf50,byid = T)
placesTransCont75 <- gIntersects(placesTransBuf75,byid = T)
placesTransCont100 <- gIntersects(placesTransBuf100,byid = T)
placesTransCont125 <- gIntersects(placesTransBuf125,byid = T)
placesTransCont150 <- gIntersects(placesTransBuf150,byid = T)
placesTransCont175 <- gIntersects(placesTransBuf175,byid = T)
placesTransCont200 <- gIntersects(placesTransBuf200,byid = T)

intersectsPlaces <- c(sum(placesTransCont25 == T),sum(placesTransCont50 == T),sum(placesTransCont75 == T),sum(placesTransCont100 == T),sum(placesTransCont125 == T),sum(placesTransCont150 == T),sum(placesTransCont175 == T),sum(placesTransCont200 == T))

plot(c(seq(25,200,25)),(intersectsPlaces-2576)/2
     ,xlab='Buffer Distance (m)'
     ,ylab='Intersections')

# merging places
# step 1 is to find the minimum population place
# step 2 is to determine if the minimum population place has a neighbor
# step 3 is to merge the minimum population location with all of its neighbors
# step 4 is to update the values for the merged polygon
# step 5 is to find the new minimum population place and repeat
# repeated until the minimum population ones can not be merged or the minimum population is above 10000

# boolean to control while loop
continueMerge <- T

# temporary spatial polygons used when merging
placesTransTemp <- placesTrans
placesTransTemp$mergeID <- NA
placesTransTemp$name2 <- placesTransTemp$NAMELSAD10

# distance for buffering
BufDist <- 50

# buffering polygons
placesTransBufTemp <- gBuffer(placesTrans # shapefile
                            ,byid=T
                            ,width=BufDist # in units of projection
)

placesTransContTemp <- gIntersects(placesTransBufTemp,byid = T)

# goal for population, will not merge areas above this population
popGoal <- 10000

# creating a sorted vector for population
popCheck <- placesTrans[['DP0010001']]

# creating list with the population sorted (smallest to largest)
# and the index of the value
popCheckSort <- as.data.frame(cbind(seq(1,length(placesTrans),1),placesTrans[['DP0010001']]))
names(popCheckSort) <- c("ind","pop")
popCheckSort$mergeID <- NA # storing the ID
popCheckSort$mergPop <- NA # storing merged population
popCheckSort$neighbors <- apply(placesTransContTemp == T,2,sum) # number of neighbors
popCheckSort$numMerges <- 0
popCheckSort$numNeighborsMerged <- NA

# creating list with the population sorted (smallest to largest)
# and the index of the value
popCheckSortTemp <- popCheckSort

# variable to control while loop
merging <- T

popCheckSortTemp2 <- popCheckSort
popCheckSortTemp2 <- popCheckSortTemp2[order(popCheckSortTemp2$pop),]


# loop for merging
while(merging){
  
  # creating ordered data frame from smallest population to largest
  popCheckSortTemp <- popCheckSortTemp[order(popCheckSortTemp$pop),]
  
  # variable controlling looping
  looping <- T
  
  # counter in the loop
  counter <- 1
  
  # looping through the data set until merged or merging not needed
  while(looping){
    
    # if population is less than the population goal
    if(popCheckSortTemp$pop[counter] < popGoal){
      
      # condition of no neighbors
      if(popCheckSortTemp$neighbors[counter] == 1){
        
        # no neighbors
        # population stays the same
        popCheckSortTemp$mergeID[counter] <- -9999
        
        # indexing counter
        counter <- counter + 1
        
        # continue to loop
        looping <- T
        print("hi-1")
      }else{
        
        # indices of neighbors
        # checking if values were already merged
        if(popCheckSortTemp$numMerges[counter] == 0){
          
          # indices that are marked as neighbors (first-order neighbors, touching)
          indsNeighbors <- which(placesTransContTemp[,popCheckSortTemp$ind[counter]] ==T)
        
          # should calculations be continued
          continueCalcs <- T
          
#         }else if(popCheckSortTemp$numMerges[counter] == 1){
#           
#           # indices that are marked as second-order neighbors (touching a neighbor)
#           indsNeighborFirst <- which(placesTransContTemp[,popCheckSortTemp$ind[counter]] ==T)
#           
#           # vector of neighbors
#           indsNeighbors <- NULL
#           
#           for(i in indsNeighborFirst){
#             
#             # finding the second order connections
#             indsNeighbors <- c(indsNeighbors,which(placesTransContTemp[,i] ==T))
#             
#           }
#           
#           # should calculations be continued
#           continueCalcs <- T
#         
        }else{
          print("hi0")
          # should calculations be continued
          # not continuing calculations beyond 2 merges
          continueCalcs <- F
          popCheckSortTemp$mergeID[counter] <- -1111
        }
        
        if(continueCalcs){
          
          # index of place
          indPlace <- popCheckSortTemp$ind[counter]
          
          # summing the population
          # finding indices in sorted matrix
          indsNeighborsSorted <- NULL
          for(i in indsNeighbors){
            indsNeighborsSorted <- c(indsNeighborsSorted,which(popCheckSortTemp$ind == i))
          }
        
          # adding the population
          popCheckSortTemp$pop[counter] <- sum(popCheckSortTemp$pop[indsNeighborsSorted],na.rm=T)
          
          # filling in NA for places that were merged
          for(i in setdiff(indsNeighbors,indPlace)){
            ind <- which(popCheckSortTemp$ind %in% i)
            popCheckSortTemp$pop[ind] <- NA
          }
          
          
          popCheckSortTemp$numMerges[counter] <- popCheckSortTemp$numMerges[counter] +1
          
          popCheckSortTemp$numNeighborsMerged[counter] <- length(indsNeighborsSorted)-1
          
          popCheckSortTemp$mergeID[indsNeighborsSorted] <- popCheckSortTemp$ind[counter] 
          print("hi1")
          # stop looping because of merges
          looping <- F
        }else{
          print("hi2")
          # continue looping
          # because of no merges
          looping <- T
          
          counter <- counter + 1
        }
      }
    }
  }

 
  print(popCheckSortTemp$pop[counter])
  #print(head(popCheckSortTemp))
  
  # condition to end the loop
  # next population has value over population goal
  if(popCheckSortTemp$pop[counter+1] >= popGoal){
    print("hi3")
    merging <- F
  }
}


placesTrans$NAMELSAD10[c(1132,969,590,959,1669,2408,7,1935,15,1172,1671,799)]








# loop over all places
for(i in 1:nrow(popCheckSort)){
  
  # condition if population is less than the population goal
  if(popCheckSort$pop[i] < popGoal){
    
    # no neighbors to merge with
    if(sum(placesTransContTemp[,popCheckSort$ind[i]]==T) == 1){
      
      popCheckSortTemp
    }
    
  }else{
    
    popCheckSortTemp$mergeID[i] <- placesTrans$GEOID10[popCheckSort$ind[i]]
    
  }
  
}



# creating a sorted vector for population
popCheck <- placesTrans[['DP0010001']]

# creating list with the population sorted (smallest to largest)
# and the index of the value
popCheckSort <- cbind(sort.list(popCheck),placesTrans[['DP0010001']][sort.list(popCheck)])

mergeIndicator <- rep(NA,nrow(popCheckSort))









# finding neighborhood based on distance
# minimum distance = d1 = 0 m
# maximum distance = d2 = 5000 m
# distance units in same units as projection
placesNeigh = dnearneigh(coordinates(placesTrans),d1=0,d2=5000)

# creating a sorted vector for population
popCheck <- placesTrans[['DP0010001']]

# creating list with the population sorted (smallest to largest)
# and the index of the value
popCheckSort <- cbind(sort.list(popCheck),placesTrans[['DP0010001']][sort.list(popCheck)])

mergeIndicator <- rep(NA,nrow(popCheckSort))

# lower limit for stopping merging
lowLimMerge <- 4000

# looping through to determine if things should be merged
for(i in 1:nrow(popCheckSort)){
  
  # -9999 if population is above value
  if(popCheckSort[i,2] > lowLimMerge){
    mergeIndicator[i] <- -9999
  }else{
    if(length(placesNeigh[[popCheckSort[i,1]]]) == 1){
      if(placesNeigh[[popCheckSort[i,1]]] == 0){
        mergeIndicator[i] <- -1111
      }else{
        mergeIndicator[i] <- popCheckSort[i,2]+placesTrans[['DP0010001']][placesNeigh[[popCheckSort[i,1]]]]
      }
    }else{
      mergeIndicator[i] <- popCheckSort[i,2]+placesTrans[['DP0010001']][placesNeigh[[popCheckSort[i,1]]]]
    }
  }
}


indsPlot <- which(mergeIndicator > 0)

hist(log10(mergeIndicator[indsPlot]))
points(log10(4000)
      ,y=0)