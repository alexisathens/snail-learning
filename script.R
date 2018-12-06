rm(list = ls())
set.seed(0)

library(raster) #geospatial
library(rgdal) #geospatial
library(dismo) #species distribution modeling, randomPoints
library(maptools) #plot spatial data
library(sp) #spatial stuff
library(mgcv) #GAMs
library(randomForest) #random forests
library(geosphere) #distm function
library(ggplot2)
library(leaps) #forward selection

#========================= Data Manipulation =============================
#import and clean data
setwd("~/Documents/UF/Fletcher Lab/SnailModel")

full <- read.csv("../NAS Giant Apple Snail Data/NAS-Data-Download.csv")
full <- full[,1:19]
#colnames(full)

#remove faulty observations
full <- full[-which(full$State=="PR"),] #remove 2 observations from Puerto Rico
full <- full[-which(full$State=="AZ"),] #remove 1 observation from Arizona
#full <- full[-which(full$State=="AZ"),] #remove x observations from Texas (read lit)

#subset into coordinates
snail <- full[,c("Longitude","Latitude")]
#colnames(snail)

#========
#download bioclim data
biodat <- raster::getData("worldclim", var="bio", res=2.5)
biodat <- dropLayer(biodat, c(3,7)) #bio3 lin dependent on bio2&7; bio7 lin dep on bio5&6

#get extent
max.lat <- ceiling(max(snail$Latitude))
min.lat <- floor(min(snail$Latitude))
max.lon <- ceiling(max(snail$Longitude))
min.lon <- floor(min(snail$Longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

#crop bioclim data to geographic extent
bioclim.dat <- crop(biodat, 1.25*geographic.extent) #to later extend sampling region

#plot raster
plot(bioclim.dat[[5]]) #plot bio6: min temp of coldest month
points(snail$Longitude, snail$Latitude, col='blue', pch=20, cex=0.75) #add observations


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
data(wrld_simpl)
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
presence <- extract(bioclim.dat, snail) #extract values from rasters
presence <- cbind(snail, presence) #include coords as predictors
absence <- extract(bioclim.dat, background)
absence <- cbind(background, absence)

#set up dataframe
dat <- data.frame()
dat <- data.frame(rbind(presence, absence))
y <- c(rep(1, nrow(presence)), rep(0, nrow(absence))) #binary response
y <- factor(y) #binary classification


#========================= Machine Learning Methods =============================

#===== classify observations into blocks
coords <- dat[,1:2] #coordinates of all points

#define many blocks
lon <- seq(from=min.lon,to=max.lon,length.out=100) #multiple of the extent
lat <- seq(from=min.lat,to=max.lat,length.out=45)
centroids <- expand.grid(lon=lon,lat=lat)
centroids$ID <- 1:nrow(centroids)

#find nearest centroid for each observation
coords$Block <- NA
for(i in 1:nrow(coords)){
  coords[i,"Block"] <- base::which.min(distm(coords[i,1:2],centroids[,1:2]))
}

#====== remove sparse blocks
#count number of points in each group
centroids$count <- NA 
for(i in 1:nrow(centroids)){
  centroids[i,"count"] <- sum(coords$Block==i)
}

#make list of centroids to be removed
rem <- which(centroids$count<=2)
centroids <- centroids[-rem,] #remove centroids with few observations
centroids$ID <- 1:nrow(centroids) #rename centroid rows
nblocks <- nrow(centroids)

#repeat blocking
for(i in 1:nrow(coords)){
  coords[i,"Block"] <- base::which.min(distm(coords[i,1:2],centroids[,1:2]))
}

#plot the classifications
gg <- ggplot(coords, aes(x=Longitude,y=Latitude,label=Block)) + 
  geom_text(check_overlap = TRUE)
gg

#count number of points in each group
centroids$count <- NA 
for(i in 1:nrow(centroids)){
  centroids[i,"count"] <- sum(coords$Block==i)
}

#stats on the blocks
nblocks #111 blocks defined
mean(centroids$count) #mean of 13
median(centroids$count) #median of 10
min(centroids$count) #min 3 observations
max(centroids$count) #max 66 observations


#====== classify observations into folds
#set.seed(0)
k <- 10 #k=10 folds

#divide blocks into k folds (from p.250 lab)
folds <- sample(1:k, nblocks, replace=TRUE)
centroids <- cbind(centroids, folds) #append corresponding fold to block number

#divide observations into folds
coords$Fold <- NA
for(i in 1:nrow(coords)){
  coords$Fold[i] <- centroids$folds[coords$Block[i]]
}
folds.obs <- coords$Fold #each observation classified into single fold


#====== forward stepwise selection
val.errors <- matrix(data=NA,nrow=ncol(dat),ncol=k)
for(i in 1:k){
  fit.fwd <- regsubsets(y[folds.obs!=i]~., data=dat[folds.obs!=i,], nvmax=ncol(dat),
                        method="forward") #get subsets for ith fold
  test.mat <- model.matrix(y[folds.obs!=i]~., data=dat[folds.obs!=i,])
  for(j in 1:ncol(dat)){
    coefi <- coef(fit.fwd,id=j)
    reg <- test.mat[,names(coefi)] %*% coefi #simple regression equation
    pred <- exp(reg)/(1+exp(reg)) #predict probabilities of group membership
    pred <- ifelse(pred>0.5,1,0) #classify to group 1 if probability > 0.5
    val.errors[j,i] <- sum(I(pred != y[folds.obs!=i]))/length(pred) #probability of misclassifying
  }
}

val.errors
#take average across columns
#then find min row -> optimal set of predictors





#======forward stepwise selection without CV
fit.fwd <- regsubsets(y~., data=dat.transformed, nvmax=ncol(dat.transformed),
                      method="forward")
fwd.sum <- summary(fit.fwd)
names(fwd.sum) #see types of decision criteria
fwd.sum$rsq # Rsq increases from 51.3% with one var to 60.7% with all vars included

#check RSS and adj-R2 stats
plot(fwd.sum$rss, xlab="Number of vars", ylab="RSS", type="l")
plot(fwd.sum$adjr2, xlab="Number of vars", ylab="Adj-Rsq", type="l")
n.fwd <- which.max(fwd.sum$adjr2) #optimal number of predictors to include is 17/19
points(n.fwd, fwd.sum$adjr2[n.fwd], col="red", cex=2, pch=20)

#do the same for Cp and BIC stats
plot(fwd.sum$cp, xlab="Number of vars", ylab="Cp", type="l")
which.min(fwd.sum$cp) #min at 16
plot(fwd.sum$bic, xlab="Number of vars", ylab="BIC", type="l")
which.min(fwd.sum$bic) #min at 8!

#plot all
plot(fit.fwd, scale="r2")

coef(fit.fwd, 8)



#===== non-spatial logistic regression
#transform predictors
dat.transformed <- base::scale(dat, center=T, scale=T)
dat.transformed <- as.data.frame(dat.transformed)

#fit logistic regression model
glm.coords <- glm(y~Longitude+Latitude, family="binomial", dat=dat.transformed) #coords only
glm.nocoords <- glm(y~., family="binomial", dat=dat.transformed[,3:19]) #all but coords
glm.few <- glm(y~bio5+bio14+bio15, family="binomial", dat=dat.transformed) #few good predictors

AIC(glm.few,glm.coords,glm.nocoords) #glm.nocoords has lowest AIC
summary(glm.nocoords)

#use train here instead

#use raster package to predict
glm.raster <- raster::predict(model=glm.nocoords, object=bioclim.dat, fun=predict)
glm.raster <- exp(glm.raster) / (1 + exp(glm.raster))
plot(glm.raster, xlab = "Longitude", ylab = "Latitude")


#GAM
#use CV for optimal # of knots
gam.nocoords <- gam(y~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+
                      bio14+bio15+bio16+bio17+bio18+bio19, 
                    family="binomial", data=dat.transformed) #all but coords


#====== random forest
#model <- y ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 +
#  bio12 + bio13 + bio14 + bio15 + bio16 + bio17
