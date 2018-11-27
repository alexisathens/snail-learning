rm(list = ls())
set.seed(0)

library(raster) #geospatial
library(rgdal) #geospatial
library(dismo) #species distribution modeling
library(maptools) #plot spatial data
library(sp) #spatial stuff
library(mgcv) #GAMs
library(randomForest) #random forests

#========================= Data Manipulation =============================
#import and clean data
setwd("~/Documents/UF/Fletcher Lab/SnailModel")

full <- read.csv("../NAS Giant Apple Snail Data/NAS-Data-Download.csv")
full <- full[,1:19]
colnames(full)

#remove faulty observations
full <- full[-which(full$State=="PR"),] #remove 2 observations from Puerto Rico
full <- full[-which(full$State=="AZ"),] #remove 1 observation from Arizona
#full <- full[-which(full$State=="AZ"),] #remove x observations from Texas (read lit)

#subset into coordinates
snail <- full[,c("Longitude","Latitude")]
colnames(snail)

#========
#download bioclim data
biodat <- getData("worldclim", var="bio", res=2.5)

#get extent
max.lat <- ceiling(max(snail$Latitude))
min.lat <- floor(min(snail$Latitude))
max.lon <- ceiling(max(snail$Longitude))
min.lon <- floor(min(snail$Longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

#crop bioclim data to geographic extent
bioclim.dat <- crop(biodat, geographic.extent)

#plot raster
plot(bioclim.dat[[6]]) #min temp of coldest month
points(snail$Longitude, snail$Latitude, col='blue', pch=20, cex=0.75)


#======
#generate pseudo-absences
#get sampling resolution from bioclim data files
bil.files <- list.files(path = "wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)

#get shape of bioclim raster
mask <- raster(bil.files[1])

#randomly sample points
background <- randomPoints(mask = mask,     #provides resolution of sampling points
                           n = nrow(snail),      #number of random points
                           ext = geographic.extent, #spatially restricts sampling
                           extf = 1.25)             #expands sampling a little bit

colnames(background) <- c("Longitude","Latitude")

#plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

#add the background points
points(background, col = "grey30", pch = 1, cex = 0.75)

#add the orignal observations
points(x = snail$Longitude, 
       y = snail$Latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

#======
#extracting values from rasters
presence <- extract(bioclim.dat, snail)
absence <- extract(bioclim.dat, background)

#set up dataframe
y <- c(rep(1, nrow(presence)), rep(0, nrow(absence))) #binary response
dat <- data.frame(cbind(y, rbind(presence, absence)))

#======
#explore colinearity
pairs(dat[,2:5], cex=0.1, fig=TRUE)



#========================= Machine Learning Methods =============================
#k-fold cv function
k <- 10
#samp <- sample(1:k, size=nrow(dat), replace=T) #but data is spatial...

#spatial clustering or deterministic lon/lat splits?
#look at quantiles of lon/lat data
coords <- rbind(snail,background)

splits <- raster::quantile(coords$Longitude,probs=seq(from=0,to=1,by=1/k))
splits <- splits[2:k] #remove 0 and 100th quantile
abline(v=splits,col="red",lty=3) 
#doesn't work because too thin of sections around high snail concentrations

#k-means spatial clustering method
k.means <- kmeans(coords, k, nstart=20)


#random forest
model <- y ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 +
  bio12 + bio13 + bio14 + bio15 + bio16 + bio17



