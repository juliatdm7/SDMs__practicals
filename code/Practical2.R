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
