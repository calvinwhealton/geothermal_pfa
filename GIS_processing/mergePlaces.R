###################
# script for merging US Census Places
# main steps are:
# 1. loading libraries
# 2. loading US Census places shapefile
# 3. determining which places are candidates for merging
# 4. mering places (sum population and area)
# 5. repeating 3 and 4 until no more merges are needed
# 6. saving merged shapefile

# at the moment the places can only be merged once
# a place can be merged with its neighboring places
# but not the neighbors of the neighbors

# loading libraries
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
library(spdep)        # other spatial processing
library(rgeos)        # for buffering places of interest and finding intersections
library(maptools)
library(xlsx)
# reading-in shapefile with US Census "places" clipped to states and region
setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization')
places <- readOGR(dsn='/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization/mergedplaces_join.shp'
                 ,layer='mergedplaces_join')

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

# distance for buffering (meters since in UTM)
BufDist <- 50

# buffering polygons
placesTransBufTemp <- gBuffer(placesTrans # shapefile
                            ,byid=T
                            ,width=BufDist # in units of projection
)

placesTransContTemp <- gIntersects(placesTransBufTemp
                                   ,byid = T)

# goal for population, will not merge places above this population to other places above this population
# will merge places above this population to places below this population
popGoal <- 10000

# creating a sorted vector for population
popCheck <- placesTrans$pop
areaLandCheck <- placesTrans[['ALAND']]
areaWaterCheck <- placesTrans[['AWATER']]
houseHoldsCheck <- placesTrans[['GeoID_1']]

# creating list with the population sorted (smallest to largest)
# and the index of the value
popCheckSort <- as.data.frame(cbind(seq(1,length(placesTrans),1),popCheck,areaLandCheck,areaWaterCheck,houseHoldsCheck))
names(popCheckSort) <- c("ind","pop","areaLand","areaWater","houseHolds")
popCheckSort$mergeID <- NA # storing the ID
#popCheckSort$mergPop <- NA # storing merged population
popCheckSort$neighbors <- apply(placesTransContTemp == T,2,sum) # number of neighbors
popCheckSort$numMerges <- 0
popCheckSort$numNeighborsMerged <- NA

# creating a temporary set of values to use in calculations
popCheckSortTemp <- popCheckSort

# variable to control while loop
merging <- T

# popCheckSortTemp2 <- popCheckSort
# popCheckSortTemp2 <- popCheckSortTemp2[order(popCheckSortTemp2$pop),]

# loop for merging
while(merging){
  
  # creating ordered data frame from smallest population to largest
  # NAs will be at the bottom of the list
  popCheckSortTemp <- popCheckSortTemp[order(popCheckSortTemp$pop),]
  
  # variable controlling looping
  looping <- T
  
  # counter in the loop
  counter <- 1
  
  # looping through the data set until merged or merging not needed
  while(looping){
    
    # checking if it has not been merged and if population is less than popGoal
    if((is.na(popCheckSortTemp$mergeID[counter])) & (popCheckSortTemp$pop[counter] < popGoal)){
      
      # checking for no neighbors
      if(popCheckSortTemp$neighbors[counter] == 1){
        
        # no neighbors
        # population and other variables stays the same
        popCheckSortTemp$mergeID[counter] <-  popCheckSortTemp$ind[counter]
        
        # indexing counter
        counter <- counter + 1
        
        # continue to loop
        looping <- T
      
      }else{ # has neighbors
        
        # indices that are marked as neighbors (first-order neighbors, touching)
        indsNeighbors <- which(placesTransContTemp[,popCheckSortTemp$ind[counter]] ==T)
        
        # index of place
        indPlace <- popCheckSortTemp$ind[counter]
        
        # index of place sorted
        indPlaceSorted <- which(popCheckSortTemp$ind == indPlace)
        
        # finding row indices in sorted dataframe
        indsNeighborsSorted <- NULL
        
        # finding the sorted data frame row indices
        for(i in indsNeighbors){
          indsNeighborsSorted <- c(indsNeighborsSorted,which(popCheckSortTemp$ind == i))
        }
        
        # removing places that were already merged from the neighbors
        indsNeighborsSortedNoNA <- indsNeighborsSorted[which(is.na(popCheckSortTemp$pop[indsNeighborsSorted]) == F)]
        
        # adding the population and other variables to the existing ones
        popCheckSortTemp$pop[counter] <- sum(popCheckSortTemp$pop[indsNeighborsSortedNoNA])
        popCheckSortTemp$areaLand[counter] <- sum(popCheckSortTemp$areaLand[indsNeighborsSortedNoNA])
        popCheckSortTemp$areaWater[counter] <- sum(popCheckSortTemp$areaWater[indsNeighborsSortedNoNA])
        popCheckSortTemp$houseHolds[counter] <- sum(popCheckSortTemp$houseHolds[indsNeighborsSortedNoNA])
        
        # assigning mergeID
        popCheckSortTemp$mergeID[indsNeighborsSortedNoNA] <- popCheckSortTemp$ind[counter]
        
        # counting the number of neighbors merged
        popCheckSortTemp$numNeighborsMerged[counter] <- length(indsNeighborsSortedNoNA)-1
        
        # filling in NA for values
        popCheckSortTemp$pop[setdiff(indsNeighborsSortedNoNA,indPlaceSorted)] <- NA
        popCheckSortTemp$areaLand[setdiff(indsNeighborsSortedNoNA,indPlaceSorted)] <- NA
        popCheckSortTemp$areaWater[setdiff(indsNeighborsSortedNoNA,indPlaceSorted)] <- NA
        popCheckSortTemp$houseHolds[setdiff(indsNeighborsSortedNoNA,indPlaceSorted)] <- NA
        
        # indexing the number of merges for the place
        popCheckSortTemp$numMerges[indsNeighborsSortedNoNA] <- popCheckSortTemp$numMerges[indsNeighborsSortedNoNA] +1
    
        # stop looping because of merges
        looping <- F
      }
    }else{
      
      # continuing to loop
      looping <- T
      
      # incrementing counter
      counter <- counter + 1
    }
  }

  print(popCheckSortTemp$pop[counter])

    # condition to end the loop
  # next population has value over population goal
  # must use counter + 1 instead of counter otherwise
  # loop will stop if the merged population is >= popGoal
  if(popCheckSortTemp$pop[counter+1] >= popGoal){
    merging <- F
  }
  print(sum(popCheckSortTemp$pop,na.rm=T))
}

# assigning the merge ID to the places that did not merge and pop > popGoal
indsHighPop <- which(popCheckSortTemp$mergeID %in% NA)
popCheckSortTemp$mergeID[indsHighPop] <- popCheckSortTemp$ind[indsHighPop]

# merging values
# assigning values to the placesTrans variable
popCheckSortTemp <- popCheckSortTemp[order(popCheckSortTemp$ind),]
placesTrans$merger <- popCheckSortTemp$mergeID
placesTrans$pop <- popCheckSortTemp$pop[popCheckSortTemp$mergeID]
placesTrans$areaLand <- popCheckSortTemp$areaLand[popCheckSortTemp$mergeID]
placesTrans$areaWater <- popCheckSortTemp$areaWater[popCheckSortTemp$mergeID]
placesTrans$houseHolds <- popCheckSortTemp$houseHolds[popCheckSortTemp$mergeID]
placesTrans$Tave_annual <- NA
placesTrans$Tmin_jan <- NA
placesTrans$Tabs_min_est <- NA
placesTrans$sfRes <- NA
placesTrans$sfAFS <- NA
placesTrans$sfASW <- NA
placesTrans$sfAER <- NA
placesTrans$sfEdu <- NA
placesTrans$sfHCSA <- NA
placesTrans$sfInf <- NA
placesTrans$sfMfg <- NA
placesTrans$sfPST <- NA
placesTrans$sfRERL <- NA
placesTrans$sfRet <- NA
placesTrans$sfWhsl <- NA
placesTrans$sfOth <- NA
placesTrans$L_network <- NA
placesTrans$GradLocal <- NA

placesTrans$Res_Bldgs <- NA
placesTrans$Res_Units <- NA
placesTrans$Com_Bldgs <- NA
placesTrans$Com_Units <- NA

placesTrans$Detached <- NA
placesTrans$Attached <- NA
placesTrans$twoto4 <- NA
placesTrans$fiveto19 <- NA
placesTrans$twentyto49 <- NA
placesTrans$fiftyplus <- NA
placesTrans$SqFtUnit <- NA
placesTrans$Tmean_jan <- NA

# assigning values to locations
for(i in 1:length(placesTrans)){
  
  ind1 <- which(placesTrans$merger %in% placesTrans$merger[i])
  ind2 <- placesTrans$GEOID[ind1]
  inds <- which(zachData$Place.Geo.Id2 %in% ind2)
  
  
  placesTrans$Tave_annual[i] <- mean(zachData$Tave_annual[inds])
  placesTrans$Tmin_jan[i] <- mean(zachData$Tmin_jan[inds])
  placesTrans$Tabs_min_est[i] <-  mean(zachData$Tabs_min_est[inds])
  placesTrans$sfRes[i] <- sum(zachData$sfRes[inds])
  placesTrans$sfAFS[i] <- sum(zachData$sfAFS[inds])
  placesTrans$sfASW[i] <- sum(zachData$sfASW[inds])
  placesTrans$sfAER[i] <- sum(zachData$sfAER[inds])
  placesTrans$sfEdu[i] <- sum(zachData$sfEdu[inds])
  placesTrans$sfHCSA[i] <- sum(zachData$sfHCSA[inds])
  placesTrans$sfInf[i] <- sum(zachData$sfInf[inds])
  placesTrans$sfMfg[i] <- sum(zachData$sfMfg[inds])
  placesTrans$sfPST[i] <- sum(zachData$sfPST[inds])
  placesTrans$sfRERL[i] <- sum(zachData$sfRERL[inds])
  placesTrans$sfRet[i] <- sum(zachData$sfRet[inds])
  placesTrans$sfWhsl[i] <- sum(zachData$sfWhsl[inds])
  placesTrans$sfOth[i] <- sum(zachData$sfOth[inds])
  placesTrans$L_network[i] <- sum(zachData$L_network[inds])
  placesTrans$GradLocal[i] <- mean(zachData$GradLocal[inds])
  
  placesTrans$Res_Bldgs[i] <- sum(zachData$Res_Bldgs[inds])
  placesTrans$Res_Units[i] <- sum(zachData$Res_Units[inds])
  placesTrans$Com_Bldgs[i] <-  sum(zachData$Com_Bldgs[inds])
  placesTrans$Com_Units[i] <- sum(zachData$Com_Units[inds])
  
  placesTrans$Detached[i] <- sum(zachData$Detached[inds])
  placesTrans$Attached[i] <- sum(zachData$Attached[inds])
  placesTrans$twoto4[i] <- sum(zachData$twoto4[inds])
  placesTrans$fiveto19[i] <- sum(zachData$fiveto19[inds])
  placesTrans$twentyto49[i] <- sum(zachData$twentyto49[inds])
  placesTrans$fiftyplus[i] <- sum(zachData$fiftyplus[inds])
  placesTrans$SqFtUnit[i] <- sum(zachData$SqFtUnit[inds])
  placesTrans$Tmean_jan[i] <- mean(zachData$Tmean_jan[inds])
  
}

writeOGR(placesTrans,dsn="placesTrans",layer="placesTrans7",driver="ESRI Shapefile")





# loading-in an excel file with other information for calculations
setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper')

# loading in the previously used data
zachData <- read.xlsx('INPUT_TABLE_STATIC_NY_PA_WV.xlsx'
                      ,header=T
                      ,sheetName='Sheet1')

length(intersect(zachData$Place.Geo.Id2,placesTrans$GEOID10))











for(i in 1:length(placesTrans)){
  
  
  
  
}


placesTransUnion <- unionSpatialPolygons(placesTrans, placesTrans$merger)




indsPlot <- which(mergeIndicator > 0)

hist(log10(mergeIndicator[indsPlot]))
points(log10(4000)
      ,y=0)