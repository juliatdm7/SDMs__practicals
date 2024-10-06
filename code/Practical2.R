# Biodiversity Under Pressure - SDMs Practicals - Practical 2                     4/10/2024

#Firstly, let's download all necessary packages and load them

install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
install.packages("predicts",dependencies=TRUE,repos="https://cloud.r-project.org")
install.packages("terra",dependencies=TRUE,repos="https://cloud.r-project.org")
library(geodata)
library(predicts)
library(terra)

#1.Downloading data from GBIF
##I've chosen to download data of the Southern marbled newt (Triturus pygmaeus), an endemism of the Iberian Peninsula that lives by the southern part of the Tajo Stream.
##Threatened by habitat loss due to desertification in Spain (mostly).
occdata <- geodata::sp_occurrence("Triturus", "pygmaeus*", geo=FALSE,removeZeros=TRUE,start=1,end=10000) #The start at the end of the name will give us the data to all other names for that species.
dim(occdata)
occdata[1:10,]

##Let's roughly plot our species' distribution data

wrld <- world(path=".")
##this function gives us an outline of the world's political boundaries. Reminder, if ever you want to know more about an R function, you can write ?function.name, e.g., ?world
plot(wrld, xlim=c(-15,5), ylim=c(32,47), col="light yellow", border="light gray")
##add the points
points(occdata$lon, occdata$lat, col='blue', pch=20)

#2. Cleaning up occurrence data
##This occurence data may have quite a few errors, so we have to clean it up before we use it.
##In order to know more or less the coordinates of the "rectangle" that contains my data, I can either ask for a summary() of the longitude and latitude or maybe plot an histogram, like following:
hist(occdata$lon)
hist(occdata$lat)
summary(occdata$lon)
summary(occdata$lat)

clean_occdata <- occdata[intersect(which(occdata$lat>35),which(occdata$lat<42)),]
clean_occdata <- clean_occdata[intersect(which(clean_occdata$lon>-10),which(clean_occdata$lon<1)),]

plot(wrld, xlim=c(-15,5), ylim=c(32,47), col="light yellow", border="light gray")
# add the points
points(clean_occdata$lon, clean_occdata$lat, col='blue', pch=20)
dim(clean_occdata)
View(clean_occdata)


#3.Downloading Worldclim data
output_dir <- "C:/Users/s2696220/OneDrive - University of Edinburgh/Biodiversity Under Pressure/Workshops_Practicals/SDMs__practicals"
bio_glob <- worldclim_global(var="bio", res=10,path=output_dir, version="2.1")

dim(bio_glob)

#we will also clip the spatraster so it only covers the spatial extent of our study species. First its longitudes then latitudes
summary(clean_occdata$lon)
summary(clean_occdata$lat)

e <- ext(-10, 10, 35, 50) #Area of our map we want our model to consider

predictors <- crop(bio_glob, e)

names(predictors) <- substring(names(predictors),11,16)
#here we're just shortening the names of predictors by taking the 11th to 16th characters.

#We can now have a look at the global climate data. Here we’ll just look at the first 9 worldclim variables
plot(predictors,1:9)


#Reminder of what each variable refers to:
##BIO1 = Annual Mean Temperature
##BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
##BIO3 = Isothermality (BIO2/BIO7) (×100)
##BIO4 = Temperature Seasonality (standard deviation ×100)
##BIO5 = Max Temperature of Warmest Month
##BIO6 = Min Temperature of Coldest Month
##BIO7 = Temperature Annual Range (BIO5-BIO6)
##BIO8 = Mean Temperature of Wettest Quarter
##BIO9 = Mean Temperature of Driest Quarter
##BIO10 = Mean Temperature of Warmest Quarter
##BIO11 = Mean Temperature of Coldest Quarter
##BIO12 = Annual Precipitation
##BIO13 = Precipitation of Wettest Month
##BIO14 = Precipitation of Driest Month
##BIO15 = Precipitation Seasonality (Coefficient of Variation)
##BIO16 = Precipitation of Wettest Quarter
##BIO17 = Precipitation of Driest Quarter
##BIO18 = Precipitation of Warmest Quarter
##BIO19 = Precipitation of Coldest Quarter

#And here we can add our species data onto a plot of climate data for the first variable.
plot(predictors,1, main = "Triturus pygmaeus vs. Annual Mean Temperature")
points(clean_occdata$lon,clean_occdata$lat, col='purple',pch=16,cex=0.2)

#4.Creating background data
##We have to do this because in this case we only have presence data and not presence/absence data
##One approach that is widely used is to sample background data from a region - this should cover locations where the species is present and is absent.

#here I'm setting the spatial extent to be broadly consistent with that of my study species (you need to make sure it is sampling from the same extent). Remember to find out how a function works you can do ?function
bg <- spatSample(predictors,5000,"random", na.rm=TRUE, as.points=TRUE,ext=e) #To know what cell number to insert, we need to test it out a little bit and see what works best. Experience helps having some intuition about what will work best.

#Here we'll plot our background points on a map of climwin variable 1 (you could change this to any of the worldclim variables)
plot(predictors, 1)
points(bg, cex=0.1)


#5.Matching occurrence and climate data
##Our final step before we can run the model is to match the climate and occurrence data - in the previous workshop this had already been done for us. Here our occurrence data are point data and our climate data are in raster (grid) format, at 10 minute resolution.

occlatlon <- cbind(clean_occdata$lon,clean_occdata$lat)
presvals <- extract(predictors, occlatlon)
#presvals is the climate data for where the species is present
backvals <- values(bg)
#backvals is the climate data for the background data
bg_lonlat <- geom(bg)
lonlats <- rbind(occlatlon, bg_lonlat[,c("x","y")])
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))
#The first column of the dataset is a vector of 1s for presences and 0s for background data.
sdmdata <- data.frame(cbind(lonlats,pb, rbind(presvals, backvals)))
#here we combine the presence and background data into a single data frame
#In sdmdata the first two columns are climate data and third column contains either a 1 (species present) or 0 (background data). The remaining columns are the corresponding worldclim variables for these location.
#In my case, all rows up until row number 3067 are "1", out of 4067

#We can also examine how colinear (i.e. correlated) predictor variables are. Highly correlated predictor variables can give rise to statistical issues.
pairs(sdmdata[,4:7], cex=0.1) #Here we just look at the correlations between the first 4 climate variables. You could extend this to look at all 19. 
