# script for loading-in general shape files

## loading-in shape files for the states
States = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# converting to UTM 17N system
NY2 <- spTransform(NY,CRS("+init=epsg:31986"))
PA2 <- spTransform(PA,CRS("+init=epsg:31986"))
WV2 <- spTransform(WV,CRS("+init=epsg:31986"))

# removing nonconverted objects to save space
rm(NY,PA,WV,States)

## loading-in shape files for the counties
Counties = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_county_WGS84_prj", stringsAsFactors=FALSE)
NY_co = Counties[which(Counties$STATEFP == "36"),]
PA_co = Counties[which(Counties$STATEFP == "42"),]
WV_co = Counties[which(Counties$STATEFP == "54"),]

# converting to UTM 17N system
NY_co2 <- spTransform(NY_co,CRS("+init=epsg:31986"))
PA_co2 <- spTransform(PA_co,CRS("+init=epsg:31986"))
WV_co2 <- spTransform(WV_co,CRS("+init=epsg:31986"))

# removing nonconverted objects to save space
rm(NY_co,PA_co,WV_co,Counties)

## loading-in city locations based on US Census Places
cities <- readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="usCensusPlaces", stringsAsFactors=FALSE)

# converting to UTM 17N system
cities2 <- spTransform(cities,CRS("+init=epsg:31986"))

# removing nonconverted objects to save space
rm(cities)

## importing places of interest
poi <- readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="placesofinterest2", stringsAsFactors=FALSE)

# converting to UTM 17 N system
poi2 <- spTransform(poi,CRS("+init=epsg:31986"))

# buffering places of interest by 5000 and 10000 m
poi2Buf5 <- gBuffer(poi2,width=5000)
poi2Buf10 <- gBuffer(poi2,width=10000)

# removing nonconverted objects to save space
rm(poi)

# importing grid locations for raster, already in UTM 17N
# fishnet <- readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="Fishnet2_label", stringsAsFactors=FALSE)


