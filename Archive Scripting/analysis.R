#analysis and stats
library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(spdep)





#let's make it spatial
species <- read.csv("Data/species.csv")

t6.sp <- readOGR("Data/GIS/T6_J.shp")
t6.sp <- readOGR("Output/T5_J.shp")
t6.sp <- readOGR("Output/T4_J.shp")


t6.sp$Width <- as.numeric(as.character((t6.sp$Width)))
t6.sp <- merge(t6.sp, species, "Species")

#spatstat
t6.ppp <- as(t6.sp, "ppp")
summary(t6.ppp)
plot(t6.ppp, which.marks = "Species")
t6.ppp$marks$Width <- as.numeric(t6.ppp$marks$Width)
print(t6.ppp$marks$Width)
plot(t6.ppp, which.marks = "Width")
marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- left_join(marks.data, species, "Species")
marks.data <- marks.data %>% dplyr::select(sp.simp, Type, Species, Width, Height)
marks(t6.ppp) <- marks.data

#make one categorical and the other numerical
t6cat <- t6.ppp
marks.data <- as.data.frame(marks(t6cat))
marks.data <- marks.data %>% dplyr::select(sp.simp, Type)
marks(t6cat) <- marks.data

t6num <- t6.ppp
marks.data <- as.data.frame(marks(t6num))
marks.data <- marks.data %>% dplyr::select(Width, Height)
marks(t6num) <- marks.data


#EDA intensity measures

ds <- density(t6.ppp)
class(ds)
plot(ds)
plot(density(split(t6.ppp)))

windows()
#EDA distance measures

pair <- distmap(t6.ppp)
plot(pair, main = "Pairwise Distances between All Distinct Pairs of Points")
plot(t6.ppp, add = TRUE)

nn <- nndist(t6.ppp)
plot(nn, main = "Nearest Neighbour Distances")
hist(nn)
plot(t6.ppp, add = TRUE)

dis <- distmap(t6.ppp)
plot(dis, main = "Empty Space Distances")
plot(t6.ppp, add = TRUE)

X <- t6.ppp
Fc <- Fest(X, "ppp")
plot(Fc)

Fest(t6.ppp, correction = "km")

plot(Fest(t6.ppp, .))

plot(Gest(t6.ppp))

plot(Kest(t6.ppp))

plot(allstats(t6.ppp))

#K function
E <- envelope(t6.ppp, Kest, nsim = 39, rank = 1)
plot(E)
Kci <- varblock(t6.ppp, Kest, nx = 3, ny = 3)
plot(Kci, iso ~ r, shade = c("loiso", "hiiso"), main = "Confidence Interval")

#G function
set.seed(120109)
r <- seq(0, 25, by = 0.05)
env.t6 <- envelope(as(t6.sp, "ppp"), fun = Gest, r = r, nrank = 2, nsim = 99)
plot(env.t6)

##marked patten analysis
summary(X)
plot(Smooth.ppp(X))
plot(markcorr(X))
plot(Emark(X))
ppm(t6cat, ~marks)

plot(alltypes(t6cat, type = "F"))
plot(alltypes(t6cat, type = "G"))

plot(Gcross(t6cat, "lar.tri", "eri.fas"))
plot(Kcross(t6cat, "lar.tri", "eri.fas"))
windows()


la <- Kmulti(t6.ppp, t6.ppp$marks$sp.simp == "lar.tri", t6.ppp$marks$sp.simp == "amb.sal", correction = "all")
plot(la)
la <- Kmulti(t6.ppp, t6.ppp$marks$Type == "shrub", t6.ppp$marks$Type == "cactus")

par(1,1)

plot(pcf(t6.ppp))
plot(frypoints(t6.ppp))
