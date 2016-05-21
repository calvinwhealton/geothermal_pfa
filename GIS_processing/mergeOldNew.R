# merging new data into old

setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization')
places <- readOGR(dsn='/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization/mergerIhopeThisWorksZachFile.shp'
                  ,layer='mergerIhopeThisWorksZachFile')

# loaing csv file
setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper')
outputData <- read.xlsx('mergedPlacesFromZachcopy.xlsx'
                  ,header=T
                  ,sheetName='OptimalTemp')

places$SLCOH <- NA
places$opTemp <- NA
places$numPlant <- NA

for(i in 1:length(places)){
  
  geoid <- places$GEOID[i]
  
  outputInd <- which(outputData$GeoID %in% geoid)
  
  places$SLCOH[i] <- outputData$LCH...MMBTU.[outputInd]
  places$opTemp[i] <- outputData$OptimalTemp[outputInd]
  places$numPlant[i] <-outputData$NumberPlants[outputInd]
  
}

writeOGR(placesTrans,dsn="placesTrans",layer="placesTrans9",driver="ESRI Shapefile")

