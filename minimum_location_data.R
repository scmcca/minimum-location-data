library(sp)
library(raster)
library(rgdal)
library(LaplacesDemon)
library(data.table)

setwd("~/Documents/R/data/TRB2021-data")
g <- read.csv("g.csv")

## add sort order for id and date

coordinates(g) <- ~lon+lat
proj4string(g) <- CRS("+init=epsg:4326")
CRS.new <- CRS("+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") 
g <- spTransform(g, CRS.new)

list <- list()
list <- split(g, interaction(g$id,g$date))

cs <- 500

rasterList <- list()
for (i in seq_along(list)){
  r <- values(
    rasterize(
      coordinates(list[[i]]),
      (raster(ext = (extent(g)+4*cs), res = cs)),
      list[[i]]$duration_s,
      fun = sum, background = 0)
  )
  rasterList[[i]] <- r/sum(r)
}

id <- vector()
date <- vector()
for (i in seq_along(list)) {
  id[i] <- unique(list[i][[1]][[1]])
  date[i] <-unique(list[i][[1]][[2]])
}

group <- list()
for (i in seq_along(id)) {
  group[[i]] <- list(
    id[[i]],
    date[[i]],
    rasterList[[i]]
  )
}

probVec <- as.data.frame(rasterList)
colnames(probVec) <- unique(
  as.vector(
    (outer(id ,date, paste, sep = "."))
  )
)

pvGroup <- list()
for (i in seq_along(unique(id))) {
  x <-   paste(sprintf("^%d.", unique(id[i])),"\\d+", sep = "")
  pvGroup[[i]] <- probVec[ , grepl(x, names(probVec), perl = TRUE) ]
}

## cumulative for each nested list
cumulativeList <- list()
temp <- list()
for (i in seq_along(pvGroup)) {
  for (j in seq_along(pvGroup[[i]])) {
    temp[[j]] <- rowSums(pvGroup[[i]][1:j]/j)
  }
  cumulativeList[[i]] <- temp
}

KLDl <- list()
KLDtemp <- list()
pypx <- list()
for (i in seq_along(cumulativeList)) {
  for (j in 1:(length(cumulativeList[[i]])-1)) {
    KLDtemp[[j]] <- KLD(
      cumulativeList[[i]][[j]],
      cumulativeList[[i]][[j+1]]
    )
    pypx[[j]] <- KLDtemp[[j]]$sum.KLD.py.px
  }
  KLDl[[i]] <- pypx
}
names(KLDl) <- unique(id)
for (i in seq_along(KLDl)) {
  for (j in seq_along(KLDl[[i]])) {
    names(KLDl[[i]]) <- tail(unique(date), -1)
  }
}

KLD <- setDT(KLDl, keep.rownames = TRUE)[]

#visualization
ncolumns <- length(KLDl[[1]])
m <- matrix(NA, length(KLDl), ncolumns)
for (i in 1:ncolumns) {
  m[ , i] <- unlist(KLD[i])
}

md <- as.data.frame(m)
mde <- vector(mode="numeric")
mdd <- vector(mode="numeric")
mdp <- vector(mode="numeric")
mdf <- vector(mode="logical")
for (i in seq_along(md)) {
  mde <- append(mde, median(md[[i]]))
}
for (i in seq_along(md)) {
  mdd <- abs(
    na.omit(
      append(mdd, (mde[i+1]-mde[1])/mde[1])
      )
    )
}
for (i in seq_along(mdd)) {
  mdp <- na.omit(
    abs(append(mdp, (mdd[i+1]-(mdd[i])))))
}
for (i in seq_along(mdp)) {
  mdf <- append(mdf,
                all(mdp[i:length(mdp)] < 0.05)==TRUE
  )
}

#plot
boxplot(m,
        xlab = "Day",
        ylab = "Kullback-Leibler divergence",
        main = "Point of marginal information gain",
        col=NULL,
        names = seq_along(KLD[[1]])+1
        )
text(length(KLD[[1]])+0.5,(max(m)-0.05*(max(m))), sprintf("Cell size: %d m", cs), cex=1, pos=2)

abline(
  v=Position(function(x) x ==TRUE , mdf)+2,
  col="red",
  lwd=1,
  lty=5)