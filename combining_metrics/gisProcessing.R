# scripts for some GIS processing
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
library(xlsx)         # for reading in xlsx tables
library(rgeos)        # for buffering places of interest
library(RColorBrewer) # R color brewer palettes for parallel axis plot
library(pracma)       # for interpolation in tables
library(vioplot)      # for violin plots
#library(GISTools)     # adding stuff to plots
library(prettymapr)
library(maps)

########################################
# calculating road length per place

roadGrid <- readOGR('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization'
                    ,layer='checker')

#plot(roadGrid)
# calculating road density per square km (km/km^2)
# original units in m/m^2
roadGrid$Density <- (roadGrid$LENGTH*1000/roadGrid$AREA)
roadGrid$Density2 <- roadGrid$Density
roadGrid$Density2[which(roadGrid$Density > 10)] <- 10
roadGrid$Density3 <- roadGrid$Density2/1000
roadGrid$lenAdj <- NA

# looping to calculate road length in the place
for(i in 1:length(unique(roadGrid$GEOID10))){
  
  q <- unique(roadGrid$GEOID10[i])
  
  # indices for geoid10
  # all the polygons that are a part of the same place
  inds <- which(roadGrid$GEOID10 %in% q)
  
  # area*density/area = length
  roadGrid$lenAdj[inds] <- sum((roadGrid$AREA[inds]*roadGrid$Density3[inds]))

}

# writing file
writeOGR(roadGrid
         ,'/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization'
         ,layer='checkerDone.shp'
         ,driver='ESRI Shapefile')

########################################
# assigning another ID to the census places

subs <- readOGR('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper'
                    ,layer='Untitled')

subs$ID2 <- subs$GEOID10

writeOGR(subs
         ,'/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization'
         ,layer='CouSub_2010Census_newID.shp'
         ,driver='ESRI Shapefile')

########################################
# reading-in data files for 2012 economic census and housing data

places <- readOGR('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization'
                    ,layer='maybeWork2')

places$no_util <- NA
places$no_mfg <- NA
places$no_rtl <- NA
places$no_trans <- NA
places$no_info <- NA
places$no_fin <- NA
places$no_real <- NA
places$no_sci <- NA
places$no_admin <- NA
places$no_educ <- NA
places$no_health <- NA
places$no_arts <- NA
places$no_food <- NA
places$no_other <- NA

setwd('/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/economicCensus')

# utilities
t1 <- read.table('EC1222A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
                 )

t_util <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]


# manufacturing
t1 <- read.table('EC1231A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_mfg <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# retail
t1 <- read.table('EC1244A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_rtl <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# transportation and warehousing
t1 <- read.table('EC1248A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_trans <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# information
t1 <- read.table('EC1251A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_info <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# finanace and insurance
t1 <- read.table('EC1252A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_fin <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# real estate and leasing
t1 <- read.table('EC1253A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_real <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# professional scientific technical servies
t1 <- read.table('EC1254A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_sci <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# administratative support waste management
t1 <- read.table('EC1256A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_admin <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# education
t1 <- read.table('EC1261A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_educ <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# healthcare
t1 <- read.table('EC1262A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_health <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# arts and entertainment
t1 <- read.table('EC1271A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_arts <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# accomodation and food
t1 <- read.table('EC1272A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_food <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]

# other
t1 <- read.table('EC1281A1.dat'
                 ,sep='|'
                 ,quote=""
                 ,header=T
)

t_other <- t1[indUse <- which(t1$ST %in% c(36,42,54)),]


# adding values to the places shapefile
# using GEOID matching between the two
# access the first value in the list because it is the coarsest category
for(i in 1:length(unique(places$GEOID10_2))){
  
  inds <- which(places$GEOID10 %in% unique(places$GEOID10)[i])
  
  # utilities
  # making matrix for the places IDs
  q <- strsplit(as.character(t_util$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% as.numeric(unique(as.character(places$GEOID10))[i]))
  #print(t2)
  places$no_util[inds] <- t_util$ESTAB[t2[1]]
  
  # manufacturing
  # making matrix for the places IDs
  q <- strsplit(as.character(t_mfg$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_mfg[inds] <- t_mfg$ESTAB[t2[1]]
  
  # retail
  # making matrix for the places IDs
  q <- strsplit(as.character(t_rtl$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_rtl[inds] <- t_rtl$ESTAB[t2[1]]
  
  # transportation
  # making matrix for the places IDs
  q <- strsplit(as.character(t_trans$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_trans[inds] <- t_trans$ESTAB[t2[1]]
  
  # information
  # making matrix for the places IDs
  q <- strsplit(as.character(t_info$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_info[inds] <- t_info$ESTAB[t2[1]]
  
  # finance
  # making matrix for the places IDs
  q <- strsplit(as.character(t_fin$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_fin[inds] <- t_fin$ESTAB[t2[1]]
  
  # real estate
  # making matrix for the places IDs
  q <- strsplit(as.character(t_real$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_real[inds] <- t_real$ESTAB[t2[1]]
  
  # science
  # making matrix for the places IDs
  q <- strsplit(as.character(t_sci$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_sci[inds] <- t_sci$ESTAB[t2[1]]
  
  # administration
  # making matrix for the places IDs
  q <- strsplit(as.character(t_admin$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_admin[inds] <- t_admin$ESTAB[t2[1]]
  
  # education
  # making matrix for the places IDs
  q <- strsplit(as.character(t_educ$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_educ[inds] <- t_educ$ESTAB[t2[1]]
  
  # health
  # making matrix for the places IDs
  q <- strsplit(as.character(t_health$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_health[inds] <- t_health$ESTAB[t2[1]]
  
  # arts
  # making matrix for the places IDs
  q <- strsplit(as.character(t_arts$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_arts[inds] <- t_arts$ESTAB[t2[1]]
  
  # food and accomodation
  # making matrix for the places IDs
  q <- strsplit(as.character(t_food$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_food[inds] <- t_food$ESTAB[t2[1]]

  # other
  # making matrix for the places IDs
  q <- strsplit(as.character(t_other$GEO_ID),"US")
  q2 <- as.data.frame(matrix(unlist(q),nrow=length(q),ncol=2,byrow=T))
  t2 <- which(as.numeric(q2[,2]) %in% unique(as.character(places$GEOID10))[i])
  places$no_other[i] <- t_other$ESTAB[t2[1]]
  
}

