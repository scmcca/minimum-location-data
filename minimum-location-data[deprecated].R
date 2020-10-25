# Computing Kullback-Leibler Divergence for Successive Days in R

# Load in appropriate libraries and set working directory
library(sp)
library(raster)
library(LaplacesDemon)
library(rgdal)
setwd("~/Documents/R/data/glh")

#Projection and Date Split
# (1) Read data for single participant and project to Quebec Lambert
pts <- read.csv("g0.csv")
coordinates(pts) <- ~lon+lat
proj4string(pts) <- CRS("+init=epsg:4326")
CRS.new <- CRS("+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

ptsn <- spTransform(pts, CRS.new)

# (2) Read data for all participants
g <- read.csv("g.csv")
coordinates(t) <- ~lon+lat
proj4string(t) <- CRS("+init=epsg:4326")
CRS.new <- CRS("+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
tn <- spTransform(t, CRS.new)

# (3) Split single participant data by date
ds <- split(ptsn, ptsn$date)

#Date-split to Probability Density Vectors
# (4) Set resolution of base raster (meters)
cs <- 100

# (5) Rasterize date-split duration values and store in list
rvl <- list()
for (i in 1:length(ds)){
  rvl[[i]] <-
    values(
      rasterize(
        coordinates(ds[[i]]),
        (raster(ext = (extent(tn)+2*cs), res = cs)),
        ds[[i]]$duration_s,
        fun = sum, background = 0)
    )
}

# (6) Convert values to probability density vectors
for (j in 1:length(rvl)){
  rvl[[j]] <-
    rvl[[j]]/sum(rvl[[j]])
}
rvl <- as.data.frame(rvl)

# (7) Produce cumulative probability density vectors
rs <- list()
for (n in 1:length(rvl)){
  rs[[n]] <- rowSums(rvl[1:(n)])/(n)
}

# (8) Calculate Kullback-Leibler Divergence for cumulative probability density
vectors
KLDl <- list()
pypx <- list()
for (d in 1:(length(rs)-1)){
  KLDl[[d]] <- KLD(rs[[d]], rs[[d+1]])
  pypx[d] <- (KLDl[[d]]$sum.KLD.py.px)
}

# (9) Write results to CSV
write.csv(unlist(pypx), file = "v100g0.csv")