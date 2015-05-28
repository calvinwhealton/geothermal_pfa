# code to analyze drilling fluid based on depth and year
# code will be stored on GitHub

## set working directory----
setwd("/Users/calvinwhealton/GitHub/geothermal/drillingfluid")

## libraries----
library(RColorBrewer)
library(scales)
library(mlogit)
## importing data----
pa <- read.csv('padatabase.csv', header=TRUE, sep=',')
wh <- read.csv('whealtonModified.csv', header=TRUE, sep=',')

### AASG PA DF----
## finding unique drilling fluids----
dfs <- unique(pa$DrillingFluid)
dfs_cnt <- rep(0,length(dfs))

for(i in 1:length(dfs)){
  dfs_cnt[i] <- length(which(pa$DrillingFluid %in% dfs[i]))
}
dfs_pct <- 100*dfs_cnt/nrow(pa)

par(mar=c(4.5,9,2,2))
barplot(dfs_pct[which(dfs_pct > 1)]
        , horiz=TRUE
        , xlab = "Percentage of Wells"
        , xlim = c(0,16))
lines(c(2.5,2.5)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(5,5)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(7.5,7.5)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(10,10)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(12.5,12.5)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(15,15)
      , c(0,20)
      , col="white"
      , lwd = 2
      , lty = 1)
text(c(rep(0,length(which(dfs_pct > 1))))
      , (c(seq(1,length(which(dfs_pct > 1)),1))-.5)*1.2
      , labels = dfs[which(dfs_pct > 1)]
      , pos =2
      , xpd = TRUE)

## correcting typos/smallest groups----
# not recorded----
norec <- c(NA, "N/A")
# empty----
empty <- c("Empty", "None", "empty"
           , "none", "Empty ", "Emprty")
# air-drilled----
ad <- c("Air Drilled", "Nitrogen And Air Drilled", "Air Dirlled"
        , "Air Drilled ", "air drilled", "Air drilled"
        , "Air Drillled")

ad_ff <- c("Formation Fluid/Air Drilled", "Air-Drilled/ Formation Fluids", "Formation Fluids/ Air Drilled"
           , "Air Drilled/ Formation Fluids ", "Air Drilled/ Formation Fluids")

ad_dry <- c("dry (air drilled)", "Dry-Air Drilled")

ad_water <- c("Air Drilled/ Water Filled", "Partially Water Filled. Air Drilled")
  
ad_damp <- c("Air Drilled(Damp)", "Damp/Air Drilled", "damp(air drilled)")

ad_other <- c("Soaped and Air Drilled","Air Drilled/Salt Gel", "air/foam")
# dry----
dry <- c("dry", "Dry", "Dru"
         , "Dry ")

dry_other <- c("Gas/Dry", "Dry Air", "Dry/Air")
# dusted----
dusted <- c("Dusted", "dusted", "Dusied")
# damp_wet----
damp <- c("Wet", "Damp", "Damp hole", "wet hole"
          , "damp", "Damp Hole")

damp_oil <- c("Damp/Oilly", "Damp/oilly", "Damp/Oily")
# soaped_foamed----
soap <- c("Soaped", "Foam", "Soap"
          , "Soap/Foam", "soap/foam", "soaped"
          , "Scaped hole", "Scaped Hole")
# air-gas----
air <- c("air",  "air ", "AIR"
         , "Air Filled", "Air Filled ", "AIr"
         , "Cable Tool/Air", "Air Filled Hole",  "Air", "Air ")

air_gas <- c("Air/Gas",   "Air/gas", "Gas/Ait"
             , "Gas/ Air", "air and gas", "air/gas"
             , "AIR/GAS","Air/ Gas", "gas/air"
             ,  "Air/Gas ", "Gas/Air","Ai/Gas"
             ,  "gas air", "air/Gas", "Gas/Air ")

gas <- c("Gas Filled", "Gas", "Gass Filled"
         , "gas", "gas ", "gas filled"
         , "gas filled", "Gas Filled ", "gas filled ")

air_gas_wd <- c("air/wet","Air Mist", "air(wet)"
                , "Damp/Air", "Damp hole, Air", "Air/Damp"
                , "Air-Damp", "air & damp", "Air and Moisture"
                , "Air,Damp", "Wet/Air", "Air & Wet"
                , "Air/wet")

air_gas_foam <- c("Air Foam", "Foam/Air", "Foam/Gas/Air"
                , "Air-Wet Foam", "Air/Soap", "Foam/air"
                , "Air/Gas/Soap", "air-foam", "air foam"
                , "Air/Foam", "Air/Gas/Soap", "Soap/Air"
                , "Air Foam ")
# air-gas-water----
air_gas_sw <- c("Sw/Air", "Air/Brine", "Gas/Salt Water"
                , "Gas/ Salt Water", "salt brine/air", "gas/air/brine"
                , "saltwater/air", "salt water/air", "air/brine"
                , "Salt Water/Air", "Air/Gas/Brine", "Air/Salt"
                , "Air/Salt Water", "Air/Salt Brine", "Salt Brine/ Air"
                , "brine and air", "brine/air", "air/salt water"
                , "Gas and Salt Water", "Salt Brine/Air", "kcl water/air")

air_gas_fw <- c("Fresh Water/Gas/Air", "Fresh Water/Gas", "Air/Fresh Water")

air_gas_w <- c("Water/Air", "Water/Air/Gas", "Gas and Water"
               , "Water/Gas", "Air/Gas/Water", "Air/Water"
               , "gas and water", "air/water", "water/gas"
               , "water/air", "Water & Air", "air & water"
               , "Gas/Water", "water-gas", "air-water"
               , "Air and Water", "H2O/Gas", "Air Water", "Mater/Gas")

air_gas_o <- c("ir/Fluid","KCL Fluid/ Air", "Fluid and Air"
               , "Fluid/Air", "Air, Fluid")

air_gas_mud <- c("Air/Mud", "Salt Mud/Air", "Air/Salt Gel"
             , "Mud/Air")

air_gas_oil <- c("Air/Oil", "Oil/Gas/Air", "Oil/Air"
                 , "Oil/Gas", "air/gas/oil", "Air/Gas/Oil"
                 , "oil-gas", "Oil/Air/Gas")
# water----
water_sw <- c("Salt Water", "Brine", "Salt brine"
        , "Saltwater", "brine", "salt water"
        , "Salt Brine", "Salt Water Brine", "Brine Water"
        , "saltwater", "salt fluid", "salty water"
        , "salt brine", "salt-brine", "brine water"
        , "brine ", "KCL/Brined Solution", "Salty Water"
        , "Salt/Brine", "Salt Water/Loaded", "Brine/Loaded"
        , "Brined", "SaltWater", "Salt Water "
        , "Salt Sat. Water", "Salt-Brine", "Salt water")

water_sw_conc <- c("KCL 2%", "KCL 5%", "3% KCL Water"
                   , "KCL Water","3% KCL", "KCL"
                   , "2% KCL Water", "4% KCL Brine", "1.5% KCL"
                   , "KCL Fluid", "Air/ KCL Water", "2%KCL Water"
                   , "Fresh/KCL", "2% KCL", "3.5% KCL Water"
                   , "calcium chloride", "Salt Sat'd"
                   , "1% KCL", "3%KCL", "CACL2"
                   , "4% KCL", "Fresh-3% KCL", "3.5% KCL"
                   , "K-CL", "3.3%KCL", "4%KCL Brine")

water_f <- c("Fresh Water", "Freshwater", "fresh water"
             , "Fresh Water", "freshwater", "Fresh Water ")

water_fs_s <- c("Foam/Saltwater", "Soap/ KCL Water", "Soap Brine"
                , "KCL/foam", "Salt Water/Foam", "Salt Water/Soap")

water_fs_o <- c("Water/Foam", "Soap/Water", "Water and Soap"
                , "foam water", "soap/water", "H2O/Foam"
                , "Soap and Water", "water/foam", "Foam/Water"
                , "Water/Soap")

water_fs_air <- c("Foam/Water/Air", "Water/Air/Foam", "Air/Foam/Water"
                  , "water/foam/air", "air/soap/water", "air,foam, water"
                  , "air/foam/water", "Water/Foam/Gas", "Air/Foam/Blackwater"
                  , "Foram/Air/Water", "Air/fFoam/ water", "Foam/Air/Water")

water_fs_air_sw <- c("air/foam/brine", "Air/Foam/Brine", "Air/Foam/Brine ")

water_o <- c("drilling water", "Water", "Waster"
             ,"water", "water level", "Water From Shallow Zone" 
             , "Water ", "water ", "water &potash"
             , "blk-water", "water&", "blackwater" 
             , "WAter", "PT Water", "Black Water"
             , "making water", "Blackwater", "pit water"
             , "Pit Water", "Wter", "water& ")

water_oil_air <- c("Water/Air/Oil", "Water/Oil/Air", "Water/Oil/Gas"
                   , "Water/oil/Air", "Water/oil/Gas", "Water/Gas/Oil"
                   , "Salt Brine/Oil/Gas", "Gas/Oil/Water", "Air/Oil/Water"
                   , "Water/Oil/gas")

water_oil <- c("Oil and Water", "Oil/Water", "Oil/Saltwater"
               ,  "oil/water", "water/oil", "Water/Oil")

water_acid <- c("Water/Acid", "acidy water", "water/acid"
                , "acid/water", "acid & water", "Water-Acid"
                , "Water & Acid", "Acid Water", "Acid/Water"
                , "water & acid", "Water and Acid", "Water/Acid"
                , "Water/Acic", "acid water", "water &acid"
                , "water & acid", "water acid", "acid-water"
                , "Acid", "Water And Acid", "Water and Acid"
                , "Aicd/Water", "Water and 7.5% HCL", "Brine and Acid"
                , "acid&water", "Water and Acid ", "water & acid "
                , "Water/Acid ", "acid &water")
# gels/polymers----
gelp_f <- c("Fresh Gel", "Fresh Gel Mud with Steel Fillings", "Water with 25 Sacks of Gel" 
        , "Fresh Water Gel",  "fresh gel", "fresh water/gel/spersene"
        , "fresh water gel", "Fresh Water Gel", "Polymer"
        , "Water Gel", "Fresh Water Gel ")

gelp_s <- c("Salt Water-Gel", "Brine/Gel", "Salt Gel Starch"
            , "Brine - Salt Gel", "Polymer/Brine", "Salt-Gel"
            , "KCL Gel", "Salt-Polymer", "4% KCL Poly Gel"
            , "Salt Polymer", "KCL+NACL+GEL", "Salt Gel"
            , "Salt Water Gel", "Salt Saturated Gel", "salt water gel"
            , "salt gel", "salt/gel", "salt gel-polymer"
            , "salt water/gel", "salt gel", "salt water & gel"
            , "salt gel/starch", "KCL/ Polymer", "KCL/Polymer"
            , "Polymer-Brine", "Salt gel", "Salt Water and Gel"
            , "Salt Gel", "KCL salt gel", "salt gel/flosal"
            , "salt polymer", "Salt Water Polymer", "Barite-Polymer-Salt"
            , "Gel/Oriskany Brine", "polymer brine", "Salt Based Polymer"
            , "KCL-Polymer", "Weighted Polymer, Calcium Chloride, Mus"
            , "salt gel ", "Salt Saturated Gel ")

gelp_z <- c("Zeogen, Emulsion", "Zeogel Impermex", "Zeogel - Dextrid", "Zeogel-Impermex"
            , "Zeogel, Impermly, Salt Base", "zeogel & flosal", "salt-zeogel"
            , "salt zeogel", "zeogel&glosal", "Zegel" 
            , "Zeogel", "zeogel & salt water", "salt water & zeogel"
            , "zeogel", "ZEDGEL/DEXTRIP", "Zeogel/Impermex"
            , "zeogel&flosal")

gelp_a <- c("Aquagel, Bariod, Caustic Soda", "Aquagel", "aqua gel"
            , "Bariod/Aquagel", "Aqua Gel", "Barold/Aquagel"
            , "Baroid/Aguagel", "Bariod/Aguagel")

gelp_o <- c("Salt Gel/Fresh Water", "gel polymer", "Poly Gel"
            , "Water-Gel", "Chemical Gel", "Gel"
            , "Gel Chem", "o-em Gel", "Gel Mud"
            , "Chem Gel", "EEO Gel", "Chem-Gel"
            , "polymer gel", "water & gel", "mylogel"
            , "XCD Polymer", "Polymer/Dextrid/Weighted", "gel")
# muds----
mud_s <- c("Salt Mud", "salt mud" , "KCL/POLYMER MUD"
        , "KCL Mud", "KCL Polymer", "Polymer Brine"
        , "Brine WBM", "Salt Mud ", "Brine Polymer"
        , "S.B. Mud", "salt mud-natural", "salty mud"
        , "mud/salt brine", "salt water-mud", "salt/mud"
        , "salt gel/foam", "salt brine/salt mud", "brine/starch"
        , "SWBM", "Salty Mud", "NaCl Mud"
        , "Brine Mud", "Bravis Dextril/KCL", "Salt Based Mud"
        , "Dextra/Bravis", "Bravis Dextrid KCL", "Hec/Strach/KCL"
        , "Bravis/ Dextro/ KCL", "KLS Free/Dextria", "Barvis/KCL"
        , "KCL?Polymer Mud", "Water/Salt Mud", "Salt/Mud"
        , "Salt mud", "Salt Polymer Mud", "Brine Polymer Mud"
        , "salt Mud", "Mud/SaltBrine", "Baravis/Brine"
        , "Baravis/KCL", "Water Base Mud", "Chem with Barite"
        , "Weighted Polymer Cacl Mud", "KCL-Gel-Soltex Mud", "Salt-Starch"
        , "Salt Starch", "Polymer Mud W/ Barite", "Salt-Polymer Mud"
        ,  "Salt Sat Mud", "Brine-Polymer Mud", "Salt Gel Mud"
        , "Cacl Mud", "SaltMud", "1.42%KCL Mud"
        , "salt polymer mud", "brine/mud", "Salt Gel ")

mud_f <- c("Fresh Mud", "Fresh Barite mud", "fresh mud"
           , "FGM", "fresh natural mud", "fresh water mud"
           , "Fresh Mud ", "Fresh Barite Mud", "Fresh Gel Mud With Steel Fillings")

mud_oil <- c("Oil", "oil based mud", "Diesel Mud Base"
             , "Oil-based Mud", "Emulsion", "Oil Based Mud"
             , "Synthetic Oil Based Mud", "Oil Emulsion", "OBM"
             , "Synthetic OBM", "Synthetic OBM", "MCA-kerosene"
             , "Synthetic Oil", "Sythetic OBM")

mud_o <- c("WBM", "Mud", "Water Based Mud"
           , "Drilling Mud", "Potash Mud", "Drlg Mud"
           , "Native Mud", "WBM", "W.B.Mud"
           , "mud ", "CBM", "Starch"
           , "Drill Mud", "starch mud", "Fresh Starch"
           , "W.B. Mud", "gel mud", "Polymer Mud"
           , "Weighted Polymer", "Bariod","Lignosulfonate"
           , "Gel Mud", "Gel and Baroid", "mud"
           , "mud-water", "native mud","mud gel"
           , "gel-mud", "water &mud", "Native Mud"
           , "Water and Mud", "Gel-Mud", "ABS-40"
           , "KemBreak", "MegaDrill", "GYP-Q-BROXEN"
           , "LBT", "WMB", "Naitve Mud"
           , "Kembreak", "Baravis/ Dextro/ KCL", "Dextra/Baravis"
           , "KLA Free/Dextria", "Baravis Dextrid KCL", "Barvis Dextril/KCL"
           , "Gel Mud ")
# formation fluids/pf----
ff_air <- c("Formation Fluid/Foam/Air", "Formation Fluid/ Air", "Formation Fluid/Air"
            , "air formation water", "formation fluid/air", "air/produced fluids"
            , "Air/Foam/Formation", "Air/Formation Fluids", "Air Formation Fluid"
            , "Formation Fluid Gas", "Air/Formation Fluid", "Air/Formation"
            , "Air/Formation Fluid", "air/formation fluid", "air&formation fluid"
            , "air-formation fluid", "formationwater/air", "air & formation fluid"
            , "Air/ Formaiton Fluid", "Formation Fluid/Aor", "Air, Formation Fluids"
            , "Air. Formation Fluid", "formation fluid&air", "air/formation fluids"
            , "Air/Fortion Fluid", "Air/Formatoin Fluid",  "Air/ Formation Fluid"
            , "Formation Fuils.Gas", "Formation Fluids/Air", "Air/Formation Water"
            , "Formation Fluid/Air/Gas", "Formation Water/Air", "Air/Gas/Formation Water"
            , "Air/Gas/Formation Fluid", "Formation Water/Gas/Air", "Formation Water/ Air/Gas"
            , "Formation Water/Gas",  "Air/Natural Formation Blackwater", "Air/Fomartion Fluid"
            , "air formation fluid", "Formation Fuilds/Air", "Air/Formation Fluid ")

ff_fs_air <- c("formation foam/air", "air-foam/formation water", "Foam/Air/Formation Fluid"
               , "Air/Foam/Formation", "Foam/Air Formation Fluids", "air,foam, and formation fluid"
               , "air, formation fluid, foam", "air, foam formation fluid",  "Formation Flui/Air/Foam"
               , "air, foam, and formation fluid", "air/foam/formation fluid", "formation fluid/air/foam"
               , "formation fluid/foam/air", "formationfluid/foam/air", "formationwater/foam/air"
               , "air/soap/formation fluid", "Formation fluid/ Air/Foam", "Foam/Air/formation Fluid"
               , "Formation Fluid/Air/Foam", "Air/Formation Fluid/Foam", "Air/Foam/Formation Fluid/Pit Fluid"
               , "Air/Foam/Formation FLuid", "Air/Soap/Formation Fluid", "Air/Foam/Formation Water"
               , "Air/Foam/Formation Fluid", "Foam/Air/Formation Fluids", "Air/Foam/Formation"
               , "Air/Foam/Formation", "Air/Foam/Formation", "air, foam, formation fluid"
               , "air,foam, and formationfluid", "Air/Foam/Formation", "Air/Foam/Formation ")

ff_fs <- c("formation fluid & foam", "formation fluid, foam", "formation water/soap"
           , "Formation Water/Soap", "Formation Fluid/Foam", "Foam/Formation Water")


ff <- c("Formation Fluids", "Making Water", "Formation Fluid"
        ,  "formation water", "formation fluid", "produced fludis"
        , "Form. Water", "Formation Water", "Formation"
        , "formation", "formation", "Fomartion Fluid"
        , "Formation fluids", "Formation", "Form. Fluid"
        , "pit/formation fluid", "formation fluid/pit fluid", "Formation Fluids "
        , "Formation fluid", "Formation Brine", "Natural Formation Blackwater"
        , "Natural Formation Fluid", "formation-water", "Formation "
        , "formation ")


ff_water_sw <- c("salt water/formation fluid", "Salt Brine/Formation Fluid")

ff_air_water <- c("air/water/formation fluid", "air,formation fluid&gel", "air/brine/formation fluid"
                  , "Foam/Air/formation FLuid/Salt Brine", "Air/Water/Formation Fluid", "Formation water ")

ff_water <- c("formation fluid/water", "water/formation fluid", "formation blackwater"
              ,"Formation Fluid/Fresh Water", "Formation water", "natural Formation Black Water/ Fresh Water"
              , "Water/Formation Fluid", "Water/Formation  Fluid", "formation fluid/water"
              , "formation fluid/water", "Water/Formation Fluid")

pf <- c("Produced Fluid", "produced Fluids", " Brine/ Produced Fluids"
        , "Produced Fluids", "Production Fluids", "produced fluids "
        , "produced fluids")

pf_air_gas <- c("air/gas/production fluids", "Air/Gas/Produced Fluids", "Air/Gas/Produced Fluid"
                , "Poduced Fluids/ Air", "Produced Fluids/Air.Gas", "Produed Fluids/Air/Gas"
                , "Produced Fluids/Gas", "Produced Flluids?air/Gas", "Produed Fluids/Air"
                , "Produced Fluid/Gas/Air", "Produced Fulids/Gas/Air", "Poduced Fluids/Air/Gas"
                , "Air/Gas/Produced Fluids", "Gas/Air/Produced Fluids", "Produced Fluid/Gas"
                , "Production Fluids/Gas/Air", "Produced/Air", "Produced Fluids/Air"
                , "Produced FLuid/Air/Gas", "Produced Fluids/Gas/Air", "Produced Fluids/Air/Gas"
                , "Production Fluids/Air/Gas", "Production Fluids/ Air/ Gas", "Production Fluids/ Air"
                , "Production Fluids/Air.Gas", "Air/Produced Fluids", "Produced Fluids/ Air "
                , "Air/Gas/Produced Fluids ", "Produced Fluid/Air/Gas")
# other fluids----
fluid <- c("Fluid",  "fluid", "Pit Fluid"
           , "Well Fluids", "Well Fluid", "fluid"
           , "Well Fluid ", "fresh fluid", "Salty Fluid"
           , "Fluid/Oil", "Fluid Water", "fluid "
           , "nautral Formation Black Water/ Fresh Water", "Foam/Air/formation Fluid/Salt Brine")

fluid_o <- c("low salt solids", "Salty", "Chem"
             , "salt chemical", "salt", "DISP.", "Oil/Bottom Well"
             , "Low Solids", "Water Base", "Water/Mud/Foam"
             , "Fresh", "Salt")

fluid_air <- c("air & fluid", "air/gas/fluid", "air, fluid"
               , "fluid/air", "air/fluid", "Air/Fluid"
               , "Air, FLuid")

fluid_air_fs <- c("Air/Fluid/Foam", "Air/Foam/Fluid", "Air/Foam/pm/Fluid"
                  , "Foam/Fluid/Air", "air/foam/fluid", "air,foam,fluid")
# other----
other <- c("Oily", "Flowline", "Fresh Air"
           , "Oli")

# CONDENSING 1
# making list from fluid to code----
df_list <- list("norec" = norec
              ,"empty" = empty
              ,"ad" = ad, "ad_ff" = ad_ff, "ad_dry" = ad_dry
                    , "ad_water"=ad_water, "ad_damp"=ad_damp, "ad_other"=ad_other
              , "dry"=dry, "dry_other"=dry_other
              , "dusted"=dusted
              , "damp"=damp, "damp_oil"=damp_oil
              , "soap"=soap
              , "air"=air, "air_gas"=air_gas, "gas"=gas
                    , "air_gas_wd"=air_gas_wd, "air_gas_foam"=air_gas_foam
              , "air_gas_fw"=air_gas_fw, "air_gas_w"=air_gas_w, "air_gas_o"=air_gas_o
                  , "air_gas_mud"=air_gas_mud, "air_gas_oil"=air_gas_oil, "air_gas_sw"=air_gas_sw
              , "water_sw"=water_sw, "water_sw_conc"=water_sw_conc, "water_f"=water_f
                    , "water_fs_s"=water_fs_s, "water_fs_o"=water_fs_o, "water_fs_air"=water_fs_air
                    , "water_fs_air_sw"=water_fs_air_sw, "water_o"=water_o, "water_oil"=water_oil
                    , "water_acid"=water_acid, "water_oil_air" =water_oil_air
              , "gelp_f"=gelp_f, "gelp_s"=gelp_s, "gelp_z"=gelp_z
                    , "gelp_a"=gelp_a, "gelp_o"=gelp_o
              , "mud_s"=mud_s, "mud_f"=mud_f, "mud_o"=mud_o
                    , "mud_oil"=mud_oil
              , "ff_air"=ff_air, "ff_fs_air"=ff_fs_air, "ff_fs"=ff_fs
                    , "ff"=ff, "ff_water_sw"=ff_water_sw, "ff_air_water"=ff_air_water
                    , "pf"=pf, "pf_air_gas"=pf_air_gas, "ff_water" = ff_water
              , "fluid"=fluid, "fluid_o"=fluid_o, "fluid_air"=fluid_air
                    , "fluid_air_fs"=fluid_air_fs
              , "other"=other)

pa$gen_name1 <- NA

for(i in 1:length(names(df_list))){
  
  names_group <- df_list[[i]]
  
  for(j in 1:as.numeric(summary(df_list)[[i]])){
    pa$gen_name1[which(pa$DrillingFluid %in% names_group[j])] <- names(df_list)[i]
  }
}

checkdefs <- pa$DrillingFluid[which(pa$gen_name1 %in% NA)]
# making barplot----
dfs_tf_names <- unique(pa$gen_name) 
dfs_tf_cnt <- rep(0,length(dfs_tf_names))
for(i in 1:length(dfs_tf_names)){
  dfs_tf_cnt[i] <- length(which(pa$gen_name %in% dfs_tf_names[i]))
}
dfs_tf_pct <- 100*dfs_tf_cnt/sum(dfs_tf_cnt)
par(mar=c(4.5,5,2,2))
barplot(dfs_tf_pct
        , horiz=TRUE
        , xlab = "Percentage of Wells"
        , xlim = c(0,20))
lines(c(2.5,2.5)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(5,5)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(7.5,7.5)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(10,10)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(12.5,12.5)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(15,15)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(17.5,17.5)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
text(c(rep(0,length(dfs_tf_names)))
     , (c(seq(1,length(dfs_tf_names),1))-.5)*1.2
     , labels = dfs_tf_names
     , pos =2
     , xpd = TRUE
     , xlim = c(0,16)
     , cex=0.5)

# CONDENSING 2
## creating mega-groups----
pa$gen_name2 <- NA
notrec <- c(norec)
noflu <- c(empty)
adall <- c(ad, ad_ff, ad_dry, ad_water, ad_damp, ad_other, dry, dusted, dry_other, damp, damp_oil)
agfs <- c(air, gas, air_gas, air_gas_wd, air_gas_foam, soap)
water_ag <- c(air_gas_fw,air_gas_w, air_gas_sw, water_fs_s, water_fs_o, water_fs_air, water_fs_air_sw)
water <- c(water_sw, water_sw_conc, water_f, water_o,water_oil
           , water_acid, water_oil_air)
mud_gels <- c(gelp_f, gelp_s, gelp_z, gelp_a, gelp_o, mud_s, mud_f, mud_o, mud_oil, air_gas_mud, air_gas_oil, air_gas_o)
ffpf <- c(ff_air, ff_fs_air, ff_fs, ff, ff_water_sw, ff_air_water
          , pf, pf_air_gas, ff_water, fluid, fluid_o, fluid_air, fluid_air_fs)
other <- c(other)

df_list2 <- list("agfs"=agfs
                , "adall" = adall
                , "water"=water
                , "mud_gels"=mud_gels
                , "ffpf"=ffpf
                , "noflu" = noflu
                , "notrec" = notrec
                , "water_ag"=water_ag
                ,"other"=other)

for(i in 1:length(names(df_list2))){
  
  names_group <- df_list2[[i]]
  
  for(j in 1:as.numeric(summary(df_list2)[[i]])){
    pa$gen_name2[which(pa$DrillingFluid %in% names_group[j])] <- names(df_list2)[i]
    
  }
}

checkdefs2 <- pa$DrillingFluid[which(pa$gen_name2 %in% NA)]
# making barplot----
dfs_conden_names <- unique(pa$gen_name2) 
dfs_conden_cnt <- rep(0,length(dfs_conden_names))
for(i in 1:length(dfs_conden_names)){
  dfs_conden_cnt[i] <- length(which(pa$gen_name2 %in% dfs_conden_names[i]))
}
dfs_conden_pct <- 100*dfs_conden_cnt/sum(dfs_conden_cnt)
par(mar=c(4.5,10.5,2,2))
barplot(dfs_conden_pct
        , horiz=TRUE
        , xlab = "Percentage of Wells"
        , xlim = c(0,40)
          )
lines(c(10,10)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(20,20)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(30,30)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
text(c(rep(0,length(dfs_conden_names)))
     , (c(seq(1,length(dfs_conden_names),1))-.5)*1.2
     , labels = c("No Fluid", "Mud/Gel/Polymer", "Not Recorded", "Air/Gas/Foam/Soap", "Air-Drilled", "Formation/Produced Fluid","Water", "Water/Air/Gas", "Other")
     , pos =2
     , xpd = TRUE
     , xlim = c(0,16)
     , cex=0.7
     )

# CONDENSING 3
# making most coarse groups----
pa$gen_name3 <- NA
all_agfsad <- c(adall, agfs)
all_mwgp <- c( water, mud_gels)
all_other <- c(other, notrec,noflu,ffpf,water_ag)

df_list3 <- list("all_agfsad" = all_agfsad
                 , "all_mwgp" = all_mwgp
                 , "all_other" = all_other)

for(i in 1:length(names(df_list3))){
  
  names_group <- df_list3[[i]]
  
  for(j in 1:as.numeric(summary(df_list3)[[i]])){
    pa$gen_name3[which(pa$DrillingFluid %in% names_group[j])] <- names(df_list3)[i]
    
  }
}
checkdefs3 <- pa$DrillingFluid[which(pa$gen_name3 %in% NA)]
# making barplot----
dfs_conden_names <- unique(pa$gen_name3) 
dfs_conden_cnt <- rep(0,length(dfs_conden_names))
for(i in 1:length(dfs_conden_names)){
  dfs_conden_cnt[i] <- length(which(pa$gen_name3 %in% dfs_conden_names[i]))
}
dfs_conden_pct <- 100*dfs_conden_cnt/sum(dfs_conden_cnt)
par(mar=c(4.5,12.5,1,1))
barplot(dfs_conden_pct
        , horiz=TRUE
        , xlab = "Percentage of Wells"
        , xlim = c(0,60)
)
lines(c(10,10)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(20,20)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(30,30)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(40,40)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(50,50)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
text(c(rep(0,length(dfs_conden_names)))
     , (c(seq(1,length(dfs_conden_names),1))-.5)*1.2
     , labels = c("Other",  "Mud/Gel/Polymer \n or Water","Air/Gas/Foam/Soap \n or Air-Drilled" )
     , pos =2
     , xpd = TRUE
     , xlim = c(0,16)
     , cex=1.2
)

# making depth-year scatterplot----
plotwells <- intersect(which(pa$DrillerTotalDepth > 0),which(pa$DrillYear > 0))

fluid_interest <- union(union(which(pa$group_condensed %in% "adall"), which(pa$group_condensed %in% "agfs")),union(which(pa$group_condensed %in% "water"), which(pa$group_condensed %in% "mud_gels")))

wells1500d <- which(pa$DrillerTotalDepth > 1500/0.3048)

plotwells_topfluids <- intersect(plotwells, fluid_interest)

plotwells_topfluids1500 <- intersect(wells1500d, plotwells_topfluids)

# assigning colors
dfcols <- brewer.pal(9, "Set1")
dfpchs <- c(15, 16, 17, 18, 0, 1, 2, 3, 4)

pa$dfcols <- NA
pa$dfpchs <- NA

for(i in 1:length(names(df_list2))){
    pa$dfcols[which(pa$group_condensed %in% names(df_list2)[i])] <- dfcols[i]
    pa$dfpchs[which(pa$group_condensed %in% names(df_list2)[i])] <- dfpchs[i]
}
# all wells----
plot(NA
     , NA
     , xlab = "Driller Total Depth (m)"
     , ylab = "Year Drilled"
     , ylim = c(1890,2015)
     , xlim = c(0,7000)
)

grid(col="lightgray"
     , lwd = 2)

points(pa$DrillerTotalDepth[plotwells]*0.3048
     , pa$DrillYear[plotwells]
     , pch = pa$dfpchs[plotwells]
     , col = pa$dfcols[plotwells]
     , cex = 0.5
     )

legend("bottomright"
       , names(df_list2)
       , pch = dfpchs
       , col = dfcols)
# top fluids-----
par(mar=c(4,4,2,2))
plot(NA
     , NA
     , xlab = "Driller Total Depth (m)"
     , ylab = "Year Drilled"
     , ylim = c(1890,2015)
     , xlim = c(0,7000)
)

grid(col="lightgray"
     , lwd = 2)

points(pa$DrillerTotalDepth[plotwells_topfluids]*0.3048
       , (pa$DrillYear + (pa$DrillMonth-1)/12)[plotwells_topfluids]
       , pch = pa$dfpchs[plotwells_topfluids]
       , col = alpha(pa$dfcols[plotwells_topfluids],0.5)
       , cex = 0.3
)

legend("bottomright"
       , c("agfs", "adall", "water", "mud_gels")
       , pch = dfpchs[1:4]
       , col = dfcols[1:4])
# top fluids 1500m-----
par(mar=c(4,4,2,2))
plot(NA
     , NA
     , xlab = "Driller Total Depth (m)"
     , ylab = "Year Drilled"
     , ylim = c(1890,2015)
     , xlim = c(1400,7000)
)

grid(col="lightgray"
     , lwd = 2)

points(pa$DrillerTotalDepth[plotwells_topfluids1500]*0.3048
       , (pa$DrillYear + (pa$DrillMonth-1)/12)[plotwells_topfluids1500]
       , pch = pa$dfpchs[plotwells_topfluids1500]
       , col = alpha(pa$dfcols[plotwells_topfluids1500],0.5)
       , cex = 0.4
)

legend("bottomright"
       , c("agfs", "adall", "water", "mud_gels")
       , pch = dfpchs[1:4]
       , col = dfcols[1:4])
# plot for spatial variation in drilling fluid----
plot(pa$LongDegree[plotwells_topfluids]
       , pa$LatDegree[plotwells_topfluids]
       , pch = pa$dfpchs[plotwells_topfluids]
       , col = alpha(pa$dfcols[plotwells_topfluids],0.5)
       , cex = 0.3
     , xlab = "Longitude"
    , ylab = "Latitude"
)

legend("bottomright"
       , c("agfs", "adall", "water", "mud_gels")
       , pch = dfpchs[1:4]
       , col = dfcols[1:4]
       )
# making plots based on decade logged----
yrs <- seq(1950,2010,by=10)
depths <- seq(500,7000,by=250)

x <- matrix(0,length(yrs),length(names(df_list2)))

for(i in 1:length(yrs)){
  
  yrins <- intersect(which(pa$LogYear >= yrs[i]), which(pa$LogYear <= yrs[i]+9))
  
  for(j in 1:length(names(df_list2))){

      x[i,j] <- length(which(pa$group_condensed[yrins] %in% names(df_list2)[[j]]))

  }
}

xrn <- matrix(0,length(yrs),length(names(df_list2)))

# renormalizing x
for(i in 1:nrow(x)){
  xrn[i,]<- x[i,]/sum(x[i,])
  
}

par(mar=(c(4,4,2,8)+0.1))
barplot(t(xrn)
        , col=dfcols
        , horiz=TRUE
        , xlab = "Proportion of Wells"
        , ylab = "Decade Logged")
par(xpd = TRUE)
legend(c(2,6)
       , names(df_list2)
       , col=dfcols
       , pch = 15
       , border="white"
       , fill = "white"
       , bty = "n")

text(x=rep(0,7)
     , y=seq(0,6,1)*1.2+.5
     , pos=2
     , labels=c('50s', '60s', '70s', '80s', '90s', '00s', '10s')
)
text(x=0.5
     , y=seq(0,6,1)*1.2+.6
     , labels = c(sum(x[1,]), sum(x[2,]), sum(x[3,]), sum(x[4,]), sum(x[5,]), sum(x[6,]), sum(x[7,]))
     , col="white")
par(xpd=FALSE)
## making plot based on depth drilled-----
y <- matrix(0,length(depths),length(names(df_list2)))

for(i in 1:length(depths)){
  
  yrins <- intersect(which(pa$DrillerTotalDepth*0.3048 >= depths[i]), which(pa$DrillerTotalDepth*0.3048 < depths[i]+250))
  
  for(j in 1:length(names(df_list2))){
    
    y[i,j] <- length(which(pa$group_condensed[yrins] %in% names(df_list2)[[j]]))
    
  }
}

yrn <- matrix(0,length(depths),length(names(df_list2)))

# renormalizing x
for(i in 1:nrow(y)){
  yrn[i,]<- y[i,]/sum(y[i,])
  
}

par(mar=(c(4,4,2,8)+0.1))
barplot(t(yrn)
        , col=dfcols
        , horiz=TRUE
        , xlab = "Proportion of Wells"
        , ylab = "Driller Total Depth")
par(xpd = TRUE)
legend(c(2,20)
       , names(df_list2)
       , col=dfcols
       , pch = 15
       , border="white"
       , fill = "white"
       , bty = "n")

text(x=rep(0,7)
     , y=seq(0,24,1)*1.2+.6
     , pos=2
     , labels=c('500-750', '750-1000', '1000-1250', '1250-1500', '1500-1750', '1750-2000','200-2250', '2250-2500', '2500-2750', '2750-3000', '3000-3250', '3250-3500', '3500-3750', '3750-4000', '4000-4250', '4250-4500', '4500-4750', '4750-5000', '5000-5250', '5250-5500', '5500-5750', '5750-6000', '6000-6250', '6250-6500', '6500-6750')
      , cex=0.5
     )
text(x=0.5
     , y=seq(0,26,1)*1.2+.7
     , labels = c(rowSums(y))
     , col="white"
     , cex=0.8)
par(xpd=FALSE)

write.csv(pa, file="padatadrillingfluid.csv", fileEncoding="utf8")

### Whealton NY and PA DF----
# CONDENSING 1
# making groups
unique(wh$FLUID_TYPE)
wh_agfs <- c("A", "AG", "AGF", "AH", "A/D"
            , "D", "SP", "VD")
wh_ad <- c("AD", "DU", "ADSW")
wh_water <- c("AW", "BR", "BST", "FW", "KCF"
              , "KCL", "KCW", "SW", "W", "WA"
              , "A&SW", "A&W")
wh_mudgel <- c("B", "BD", "BG", "BP", "BZ"
               , "C", "CG", "DX", "EM", "FC"
               , "FG", "FGM", "FM", "FP", "G"
               , "GB", "GL", "GM", "HTC", "KCG"
               , "KCP", "LS", "LSM", "LSS", "M"
               , "NV", "O", "OE", "P", "PM"
               , "QG", "SBM", "SC", "SG", "SM"
               , "SPM", "ST", "WBG", "WBM", "WP"
               , "XP", "ZG", "KCM", "KC")
wh_ffpf <- c("AFF", "FF", "AG&FF")
wh_fluid <- c("F", "FL", "OGW", "S")
wh_other <- c("BL", NA)
wh_nf <- c("E", "N")

wh$gen_name1 <- NA
df_listwh <- list("agfs" = wh_agfs
                  ,"adall" = wh_ad
                  , "water"=wh_water
                  , "mud_gels"=wh_mudgel
                  , "ffpf"=wh_ffpf
                  , "noflu" = wh_nf
                  , "other"=c(wh_other, wh_fluid)
                )

for(i in 1:length(names(df_listwh))){
  
  names_group <- df_listwh[[i]]
  
  for(j in 1:as.numeric(summary(df_listwh)[[i]])){
    wh$gen_name1[which(wh$FLUID_TYPE %in% names_group[j])] <- names(df_listwh)[i]
  }
}

checkdefs <- wh$FLUID_TYPE[which(wh$gen_name1 %in% NA)]
# making barplot----
dfs_conden_names <- unique(wh$gen_name1) 
dfs_conden_cnt <- rep(0,length(dfs_conden_names))
for(i in 1:length(dfs_conden_names)){
  dfs_conden_cnt[i] <- length(which(wh$gen_name1 %in% dfs_conden_names[i]))
}
dfs_conden_pct <- 100*dfs_conden_cnt/sum(dfs_conden_cnt)
par(mar=c(4.5,10.5,2,2))
barplot(dfs_conden_pct
        , horiz=TRUE
        , xlab = "Percentage of Wells"
        , xlim = c(0,50)
)
lines(c(10,10)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(20,20)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(30,30)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(40,40)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
lines(c(50,50)
      , c(0,60)
      , col="white"
      , lwd = 2
      , lty = 1)
text(c(rep(0,length(dfs_conden_names)))
     , (c(seq(1,length(dfs_conden_names),1))-.5)*1.2
     , labels = c("Water","Mud/Gel/Polymer","Other","Air-Drilled","Air/Gas/Foam/Soap",  "No Fluid","Formation/Produced Fluid" )
     , pos =2
     , xpd = TRUE
     , xlim = c(0,16)
     , cex=0.7
)

## CONDENSING 2
#
#################################
# fitting multinomial logit model
#################################

# cleaning data
pa$yr_diff <- abs(pa$DrillYear - pa$LogYear)

ind_good_yr <- union(which(pa$yr_diff %in% 0), which(pa$yr_diff %in% 1))

intersect(ind_good_yr, which(pa$DrillYear %in% 1900)) # check of 1900's dropped

pa$depth_m <- pa$DrillerTotalDepth*0.3048

pa_mat <- matrix(0,length(ind_good_yr),2)
pa_mat[,1] <- pa$depth_m[ind_good_yr]
pa_mat[,2] <- pa$LogYear[ind_good_yr]
pa_clean <- data.frame(pa_mat)
colnames(pa_clean) <- c("depth", "logYear")
pa_clean$cat <- pa$coarse_gp[ind_good_yr]

pa_ml <- mlogit.data(pa_clean, choice = "cat")


