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
library(ggplot2) #plot graphs
library(leaps) #forward selection
#library(caret) #train function
#library(MASS) #stepAIC function
library(splines) #ns function in GAM
library(glmnet) #lasso
library(pROC) #ROC curves

library(jagsUI) #JAGS
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
set.seed(0)
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
  geom_text(check_overlap = TRUE) + ggtitle("Spatially blocked observations")
gg

#count number of points in each group
centroids$count <- NA 
for(i in 1:nrow(centroids)){
  centroids[i,"count"] <- sum(coords$Block==i)
}

#stats on the blocks
nblocks #116 blocks defined
mean(centroids$count) #mean of 12.5
median(centroids$count) #median of 9
min(centroids$count) #min 3 observations
max(centroids$count) #max 48 observations


#====== classify observations into folds
set.seed(0)
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


#===== simple logistic regression
set.seed(0)
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  glm.fit <- glm(y[folds.obs!=i]~., family="binomial", dat=dat[folds.obs!=i,3:19]) #all but coords
  pred <- predict(glm.fit, newdata=dat[folds.obs==i,3:19], type="response") #group membership probs
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
glm.misclass <- mean(val.errors)
glm.misclass #misclassified observations 12% of time

glm.roc <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(glm.roc)
glm.roc$auc #0.8713


#===== forward feature selection
set.seed(0)
val.errors <- matrix(data=NA,nrow=ncol(dat),ncol=k)

for(i in 1:k){
  fit.fwd <- regsubsets(y[folds.obs!=i]~., data=dat[folds.obs!=i,], nvmax=ncol(dat),
                        method="forward") #order best predictors for ith fold
  for(j in 1:ncol(dat)){
    coefi <- coef(fit.fwd, id=j) #get names of best j predictors
    best.pred <- names(coefi)[-1] #get rid of intercept
    #write new formula
    f <- as.formula(paste("y[folds.obs!=i] ~", paste(best.pred, collapse = " + ")))
    #use glm to predict class
    glm.fit <- glm(f, family="binomial", dat=dat[folds.obs!=i,])
    pred <- predict(glm.fit, newdata=dat[folds.obs==i,], type="response") #get probs
    pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
    val.errors[j,i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  }
}
avg.val.errors <- apply(val.errors,1,mean) #take average performance across k folds (cols)
avg.val.errors

opt.pred <- base::which.min(avg.val.errors) #select set of predictors with lowest val error
opt.pred #optimal to use 14/19 predictors
avg.val.errors[opt.pred] #misclassification rate of 11.11%


#what are these predictors? which are the best?

fit.fwd <- regsubsets(y~., data=dat, nvmax=ncol(dat), method="forward")
fwd.sum <- summary(fit.fwd)
fwd.sum$rsq # Rsq increases from 48.9% with one var to 58.7% with 17 best/all vars included
#plot rsq
plot(fwd.sum$rsq, xlab="Number of vars", ylab="Rsq", type="l")
points(opt.pred, fwd.sum$rsq[opt.pred], col="red", cex=2, pch=20)


#===== GAM
set.seed(0)

#==natural splines for best predictor
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  #fit natural splines for 5 best predictors
  gam.fit.1 <- gam(y[folds.obs!=i]~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+
                   bio13+bio14+bio15+bio16+bio17+ns(bio18,5)+bio19, family="binomial", 
                 data=dat[folds.obs!=i,])
  pred <- predict(gam.fit.1, newdata=dat[folds.obs==i,3:19], type="response") #group membership probs
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
gam.misclass.1 <- mean(val.errors) 
gam.misclass.1 #misclassified observations 12% of time

gam.roc.1 <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(gam.roc.1)
gam.roc.1$auc #0.8812


#==natural splines for best 3 predictors
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  #fit natural splines for 3 best predictors
  gam.fit.3 <- gam(y[folds.obs!=i]~bio1+bio2+bio4+bio5+bio6+ns(bio8,5)+bio9+bio10+bio11+bio12+
                     ns(bio13,5)+bio14+bio15+bio16+bio17+ns(bio18,5)+bio19, family="binomial", 
                   data=dat[folds.obs!=i,])
  pred <- predict(gam.fit.3, newdata=dat[folds.obs==i,3:19], type="response") #group membership probs
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
gam.misclass.3 <- mean(val.errors) 
gam.misclass.3 #misclassified observations 12% of time

gam.roc.3 <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(gam.roc.3)
gam.roc.3$auc #0.8812


#==natural splines for best 5 predictors
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  #fit natural splines for 5 best predictors
  gam.fit.5 <- gam(y[folds.obs!=i]~ns(bio1,5)+ns(bio2,5)+bio4+bio5+bio6+ns(bio8,5)+bio9+bio10+
                     bio11+bio12+ns(bio13,5)+bio14+bio15+bio16+bio17+ns(bio18,5)+bio19, 
                   family="binomial", data=dat[folds.obs!=i,])
  pred <- predict(gam.fit.5, newdata=dat[folds.obs==i,3:19], type="response") #group membership probs
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
gam.misclass.5 <- mean(val.errors) 
gam.misclass.5 #misclassified observations 13% of time

gam.roc.5 <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(gam.roc.5)
gam.roc.5$auc #0.8812


#==natural splines for best 8 predictors
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  #fit natural splines for 8 best predictors
  gam.fit.8 <- gam(y[folds.obs!=i]~ns(bio1,5)+ns(bio2,5)+bio4+ns(bio5,5)+ns(bio6,5)+ns(bio8,5)+
                     bio9+bio10+ns(bio11,5)+bio12+ns(bio13,5)+bio14+bio15+bio16+bio17+ns(bio18,5)+
                     bio19, family="binomial", data=dat[folds.obs!=i,])
  pred <- predict(gam.fit.8, newdata=dat[folds.obs==i,3:19], type="response") #group membership probs
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
gam.misclass.8 <- mean(val.errors) 
gam.misclass.8 #misclassified observations 13.5% of time

gam.roc.8 <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(gam.roc.8)
gam.roc.8$auc #0.8469


#====== random forest
set.seed(0)
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()

for(i in 1:k){
  forest.fit <- randomForest(y[folds.obs!=i]~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+
                               bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19,
                             data=dat[folds.obs!=i,3:19]) #all but coords
  pred <- predict(forest.fit, newdata=dat[folds.obs==i,3:19], type="response") #predict category
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
forest.misclass <- mean(val.errors)
forest.misclass #misclassified observations 12% of time

forest.roc <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(forest.roc)
forest.roc$auc #0.8362


#====== lasso shrinkage method
set.seed(0)
val.errors <- rep(NA, k)
pred.vals <- actual.vals <- numeric()
grid=10^seq(10,-2,length=100) #lambda values for lasso

for(i in 1:k){
  lasso.fit <- glmnet(x=as.matrix(dat[folds.obs!=i,]),y=as.matrix(y[folds.obs!=i]), alpha=1,
                       lambda=grid) #calculate with many tuning parameters
  min.lambda <- min(lasso.fit$lambda)
  pred <- predict(lasso.full, s=min.lambda, newx=as.matrix(dat[folds.obs==i,]))
  pred <- ifelse(pred>0.5,1,0) #classify to presence (1) if probability > 0.5
  val.errors[i] <- sum(I(pred != y[folds.obs==i]))/length(pred) #prob of misclassifying
  
  #for roc curve
  pred.vals <- c(pred.vals, pred)
  actual.vals <- c(actual.vals, y[folds.obs==i])
}
lasso.misclass <- mean(val.errors) 
lasso.misclass #misclassified observations 12.3% of time

plot(lasso.fit) #see effect of different tuning parameters

lasso.roc <- roc(as.numeric(pred.vals), as.numeric(actual.vals))
plot(lasso.roc)
lasso.roc$auc #0.8689


#====== Bayesian binary probit
prob.fit <- glm(y~., family=binomial(link="probit"), 
                dat=dat[,3:19]) #all but coords

sum.prob <- summary(prob.fit)

#run JAGS
set.seed(0)
nobs <- length(dat)

#MCMC settings 
ni <- 10000 #number of iterations
nt <- 5 #interval to thin 
nb <- 5000#number of iterations to discard as burn-in
nc <- 3 #number of chains

#set parameters to monitor
params=c("b5","b12", "b14", "b15", "b16", "b18", "b19")

#set initial values for the y's
yinit=as.numeric(y)
yinit[]=NA
inits<-function(){
  list(y=yinit)
}

#run the model
model.j=jags(model.file="snail_JAGS.R",
            parameters.to.save=params,inits=inits,
            data=list(y=as.numeric(y), x1=dat$bio5,
                      nobs=nrow(dat)),
            n.chains=nc,n.burnin=nb,n.iter=ni,n.thin=nt,DIC=TRUE)



#====== plot ROC curves
plot(glm.roc, col = 1, lty = 2, main = "ROC")
plot(gam.roc.3, col = 2, lty = 3, add = TRUE)
plot(gam.roc.5, col = 3, lty = 3, add = TRUE)
plot(gam.roc.8, col = 4, lty = 3, add = TRUE)
plot(forest.roc, col = 5, lty = 3, add = TRUE)
plot(lasso.roc, col = 6, lty = 3, add = TRUE)
legend(0.4, 0.7, legend=c("GLM", "GAM 3 NS", "GAM 5 NS", "GAM 8 NS", "R.  Forest", "Lasso"),
       col=c(1,2,3,4,5,6), lty=c(2,3,3,3,3,3), cex=0.8, box.lty=0)


#====== use raster package to predict and plot snail distribution
glm.raster <- raster::predict(model=glm.fit, object=bioclim.dat, fun=predict)
glm.raster <- exp(glm.raster) / (1 + exp(glm.raster))
plot(glm.raster, xlab = "Longitude", ylab = "Latitude", 
     main="Projected distribution of apple snail")

#====== plot best predictor
plot(bioclim.dat[[16]], main="Precipitation of warmest quarter") #plot bio18
points(snail$Longitude, snail$Latitude, col='blue', pch=20, cex=0.75) #add observations

plot(bioclim.dat[[6]], main="Mean temp of wettest quarter") #plot bio8
points(snail$Longitude, snail$Latitude, col='blue', pch=20, cex=0.75) #add observations

plot(bioclim.dat[[11]], main="Precipitation of wettest month") #plot bio13
points(snail$Longitude, snail$Latitude, col='blue', pch=20, cex=0.75) #add observations


